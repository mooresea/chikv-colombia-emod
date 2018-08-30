import os

from COMPS.Data.AssetCollection import AssetCollection as COMPSAssetCollection
from simtools.AssetManager.AssetCollection import AssetCollection
from simtools.AssetManager.FileList import FileList
from simtools.SetupParser import SetupParser
from simtools.Utilities.COMPSUtilities import get_asset_collection
from simtools.Utilities.General import init_logging

logger = init_logging("SimulationAssets")


class SimulationAssets(object):
    """
    This class represents a set of AssetCollection objects that together define all files needed by a simulation.
    """
    class InvalidCollection(Exception): pass
    class NotPrepared(Exception): pass

    EXE = 'exe'
    DLL = 'dll'
    INPUT = 'input'
    LOCAL = 'local'
    MASTER = 'master'
    PYTHON = 'python'

    COLLECTION_TYPES = [DLL, EXE, INPUT, PYTHON]
    SETUP_MAPPING = {DLL: 'dll_root', EXE: 'exe_path', INPUT: 'input_root', PYTHON: 'python_path'}

    def __init__(self):
        self.collections = {}
        self.base_collections = {}
        self.experiment_files = FileList()
        self.prepared = False
        self.master_collection = None
        self.paths = {}

    def _get_path(self, path_type):
        return self.paths.get(path_type, None) or SetupParser.get(self.SETUP_MAPPING[path_type])

    def _set_path(self, path_type, path, is_dir=False):
        if not os.path.exists(path) or (is_dir and not os.path.isdir(path)):
            raise Exception("The path specified in {} does not exist ({})".format(self.SETUP_MAPPING[path_type], path))

        # Remove the base collection
        self.base_collections[path_type] = None

        # Set the path
        self.paths[path_type] = path

    @property
    def exe_name(self):
        return os.path.basename(self.exe_path)

    @property
    def exe_path(self):
        if self.master_collection:
            for f in self.master_collection.asset_files_to_use:
                if f.file_name.endswith("exe"):
                    return f.file_name

        # If we do have a EXE collection, the exe_path should come from there
        if self.EXE in self.base_collections and self.base_collections[self.EXE]:
            for f in self.base_collections[self.EXE].asset_files_to_use:
                if f.file_name.endswith("exe"):
                    return f.file_name
            raise Exception("The asset collection ({}) used for executable does not contain any executable..."
                            .format(str(self.base_collections[self.EXE].base_collection.id)))

        return self._get_path(self.EXE)

    @exe_path.setter
    def exe_path(self, exe_path):
        self._set_path(self.EXE, exe_path)

    @property
    def input_root(self):
        return self._get_path(self.INPUT)

    @input_root.setter
    def input_root(self, input_root):
        self._set_path(self.INPUT, input_root, is_dir=True)

    @property
    def dll_root(self):
        return self._get_path(self.DLL)

    @dll_root.setter
    def dll_root(self, dll_root):
        self._set_path(self.DLL, dll_root, is_dir=True)

    @property
    def python_path(self):
        return self._get_path(self.PYTHON)

    @python_path.setter
    def python_path(self, python_path):
        self._set_path(self.PYTHON, python_path, is_dir=True)

    def __contains__(self, item):
        for col in self.collections.values():
            if item in col: return True
        return False

    @property
    def collection_id(self):
        if not self.prepared or not self.master_collection:
            return None
        return self.master_collection.collection_id

    def set_base_collection(self, collection_type, collection):
        # Make sure we have the good inputs
        if collection_type not in self.COLLECTION_TYPES and collection_type != self.MASTER:
            raise self.InvalidCollection("Collection type %s is not supported..." % collection_type)

        if not collection:
            raise self.InvalidCollection("No collection provided in set_input_collection.")

        # If the collection given is not already an AssetCollection -> retrieve
        if not isinstance(collection, COMPSAssetCollection):
            collection_id = collection
            with SetupParser.TemporarySetup('HPC'):
                collection = get_asset_collection(collection_id)

            if not collection:
                raise self.InvalidCollection("The input collection '%s' provided could not be found on COMPS." % collection_id)

        if collection_type == self.MASTER:
            self.master_collection = AssetCollection(base_collection=collection)
            for t in self.SETUP_MAPPING: self.SETUP_MAPPING[t] = ""

        else:
            self.base_collections[collection_type] = AssetCollection(base_collection=collection)

        # For the DLL filter out the exe as we usually have exe+dll in the same collection
        if collection_type == self.DLL:
            self.base_collections[self.DLL].asset_files_to_use = [a for a in self.base_collections[self.DLL].asset_files_to_use if not a.file_name.endswith('exe')]

    def prepare(self, config_builder):
        """
        Calls prepare() on all unprepared contained AssetCollection objects.
        :return: Nothing
        """
        location = SetupParser.get("type")
        self.collections = {}
        self.create_collections(config_builder)

        for collection in self.collections.values():
            if not collection.prepared:
                collection.prepare(location=location)

        # Sort the collections to first gather the assets from base collections and finish with the locally generated
        # Use a set to remove duplicates
        sorted_collections = sorted(set(self.collections.values()), key=lambda x: x.base_collection is None)

        # Gather the collection_ids from the above collections now that they have been prepared/uploaded (as needed)
        # and generate a 'super AssetCollection' containing all file references.
        asset_files = {}
        for collection in sorted_collections:
            if location == 'LOCAL':
                for asset in collection.asset_files_to_use:
                    asset_files[(asset.file_name, asset.relative_path)] = asset
            else:
                if collection.collection_id is not None: # None means "no files in this collection"
                    # Add the assets to the asset_files
                    # Make sure we have a key (file_name, path) to make sure we override assets
                    for asset in collection.comps_collection.assets:
                        asset_files[(asset.file_name, asset.relative_path)] = asset

        # Delete collections that are None (no files)
        self.collections = {cid:collection for cid, collection in self.collections.items() if collection.collection_id}

        logger.debug("Creating master collection with %d files" % len(asset_files))
        if len(asset_files) != 0:
            self.master_collection = AssetCollection(remote_files=asset_files.values())
            self.master_collection.prepare(location=location)
        self.prepared = True

    def create_collections(self, config_builder):
        location = SetupParser.get("type")
        for collection_type in self.COLLECTION_TYPES:
            # Dont do anything if already set
            if collection_type in self.collections and self.collections[collection_type]: continue

            # Check if we are running LOCAL we should not have any base collections
            if location == "LOCAL" and collection_type in self.base_collections:
                print("The base_collection of type %s was specified but you are trying to run a LOCAL experiment.\n"
                          "Using COMPS collection with LOCAL experiments is not supported. The collection will be ignored..." % collection_type)

            if location == "HPC":
                # If we already have the master collection set -> set it as collection for every types
                # Bypass the python type for now
                if self.master_collection and collection_type != self.PYTHON:
                    self.collections[collection_type] = self.master_collection
                    continue

                if collection_type not in self.base_collections:
                    base_id = SetupParser.get('base_collection_id_%s' % collection_type, None)
                    if base_id: self.set_base_collection(collection_type, base_id)

            base_collection = self.base_collections.get(collection_type, None)
            if not base_collection:
                files = self._gather_files(config_builder, collection_type)
                if files: self.collections[collection_type] = AssetCollection(local_files=files)
            else:
                self.collections[collection_type] = base_collection

        # If there are manually added files -> add them now
        if self.experiment_files:
            self.collections[self.LOCAL] = AssetCollection(local_files=self.experiment_files)

    def _gather_files(self, config_builder, collection_type):
        """
        Identifies local files associated with the given collection_type.
        :param config_builder: A ConfigBuilder object associated with this process
        :param collection_type: one of cls.COLLECTION_TYPES
        :return: A FileList object
        """
        file_list = None
        if collection_type == self.EXE:
            exe_path = self.exe_path
            if not exe_path: return None
            file_list = FileList(root=os.path.dirname(exe_path),
                                 files_in_root=[os.path.basename(exe_path)],
                                 ignore_missing=config_builder.ignore_missing)
        elif collection_type == self.INPUT:
            # returns a Hash with some items that need filtering through
            input_files = config_builder.get_input_file_paths()
            if input_files:
                file_list = FileList(root=self.input_root,
                                     files_in_root=input_files,
                                     ignore_missing=config_builder.ignore_missing)
        elif collection_type == self.DLL:
            dll_relative_paths = config_builder.get_dll_paths_for_asset_manager()
            if dll_relative_paths:
                file_list = FileList(root=self.dll_root,
                                     files_in_root=dll_relative_paths,
                                     ignore_missing=config_builder.ignore_missing)
        elif collection_type == self.PYTHON:
            file_list = FileList(root=self.python_path, relative_path="python") if self.python_path else None
        else:
            raise Exception("Unknown asset classification: %s" % collection_type)
        return file_list

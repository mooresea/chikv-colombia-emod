import argparse
import os
from calibtool import commands_args
from simtools.SetupParser import SetupParser
import simtools.Utilities.Initialization as init
from calibtool.ResampleManager import ResampleManager



def get_calib_manager(args, unknownArgs, force_metadata=False):
    mod = args.loaded_module

    if not hasattr(mod, 'run_calib_args') or 'calib_manager' not in mod.run_calib_args:
        warning_note = \
            """
            /!\\ WARNING /!\\ Suggest to set calib_manager within run_calib_args like the following:
            
                run_calib_args = {'calib_manager': calib_manager}
            """
        print(warning_note)
        manager = mod.calib_manager
    else:
        manager = mod.run_calib_args['calib_manager']

    # Atiye: uncomment this block if quiet is not being honored down below in the code. If it isn't we should fix
    #  that anyway. This message will not self destruct, ever.
    # try:
    #     if args.quiet:
    #         mod.run_sim_args['quiet'] = True
    # except:
    #     pass

    # Update the SetupParser to match the existing experiment environment/block if force_metadata == True
    if force_metadata:
        exp = manager.get_experiment_from_iteration(iteration=args.iteration, force_metadata=force_metadata)
        SetupParser.override_block(exp.selected_block)
    return manager


def get_resamplers(args, unknownArgs, force_metadata=False):
    mod = args.loaded_module

    if not hasattr(mod, 'run_calib_args') or 'resamplers' not in mod.run_calib_args:
        warning_note = \
            """
            /!\\ WARNING /!\\ Required to set resamplers within run_calib_args like the following:

                run_calib_args = {'resamplers': [xxx_resampler1, yyy_resampler2, ...]}
            """
        print(warning_note)
        exit()
    else:
        resamplers = mod.run_calib_args['resamplers']

    # Update the SetupParser to match the existing experiment environment/block if force_metadata == True
    if force_metadata:
        exp = resamplers[0].calib_manager.get_experiment_from_iteration(iteration=args.iteration, force_metadata=force_metadata)
        SetupParser.override_block(exp.selected_block)

    return resamplers


def run(args, unknownArgs):
    manager = get_calib_manager(args, unknownArgs)
    manager.run_calibration()


def resample(args, unknownArgs):
    # step 1: Get the resamplers and calibration manager
    resamplers = get_resamplers(args, unknownArgs)
    calib_manager = get_calib_manager(args, unknownArgs)
    resample_manager = ResampleManager(steps=resamplers, calibration_manager=calib_manager,
                                       restart_at_step=args.restart)

    # step 2: Resample!
    resample_manager.resample_and_run()


def resume(args, unknownArgs):
    if args.iter_step:
        if args.iter_step not in ['commission', 'analyze', 'plot', 'next_point']:
            print("Invalid iter_step '%s', ignored." % args.iter_step)
            exit()
    manager = get_calib_manager(args, unknownArgs, force_metadata=True)
    manager.resume_calibration(args.iteration, iter_step=args.iter_step)


def cleanup(args, unknownArgs):
    manager = args.loaded_module.calib_manager
    # If no result present -> just exit
    if not os.path.exists(os.path.join(os.getcwd(), manager.name)):
        print('No calibration to delete. Exiting...')
        exit()
    manager.cleanup()


def kill(args, unknownArgs):
    manager = args.loaded_module.calib_manager
    manager.kill()


def main():
    parser = argparse.ArgumentParser(prog='calibtool')
    subparsers = parser.add_subparsers()

    # 'calibtool run' options
    commands_args.populate_run_arguments(subparsers, run)

    # 'calibtool resample' options
    commands_args.populate_resample_arguments(subparsers, resample)

    # 'calibtool resume' options
    commands_args.populate_resume_arguments(subparsers, resume)

    # 'calibtool cleanup' options
    commands_args.populate_cleanup_arguments(subparsers, cleanup)

    # 'calibtool kill' options
    commands_args.populate_kill_arguments(subparsers, kill)

    # run specified function passing in function-specific arguments
    args, unknownArgs = parser.parse_known_args()
    # This is it! This is where SetupParser gets set once and for all. Until you run 'dtk COMMAND' again, that is.
    init.initialize_SetupParser_from_args(args, unknownArgs)
    args.func(args, unknownArgs)

if __name__ == '__main__':
    main()

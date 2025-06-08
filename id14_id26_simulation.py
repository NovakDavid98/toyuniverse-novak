                print("INFO: Successfully saved time-series data to .npy files.", flush=True)
            print_memory_usage("After NPY Save Loop")

        except Exception as e_outer_block:
            print(f"ERROR: A critical error occurred in the NPY saving block (outside individual file saves): {e_outer_block}", flush=True)
            # save_error = True # This would be redundant as the block itself failed
    else:
        print("INFO: Did not proceed with saving .npy files due to missing or invalid time-series lists (initial check failed).", flush=True)
        print_memory_usage("After NPY Saving Block (skipped saves due to initial check)")

    # This comment ensures no other code is accidentally appended here by mistake.

    print("\n--- Script End ---") # Final end marker

if __name__ == '__main__':
    profiler = cProfile.Profile()
    profiler.enable()
    
    run_profiler(args) # Call the main function that now contains all the script logic
    
    profiler.disable()
    
    s = io.StringIO()
    # Sort stats by cumulative time
    ps = pstats.Stats(profiler, stream=s).sort_stats('cumulative')
    ps.print_stats(20) # Print top 20 functions by cumulative time
    
    print("\n--- cProfile Stats (Top 20 by cumulative time) ---")
    print(s.getvalue())
    
    # Save full stats to a file for more detailed analysis if needed
    profile_output_filename = os.path.join(MAIN_OUTPUT_DIR if 'MAIN_OUTPUT_DIR' in globals() and MAIN_OUTPUT_DIR else '.', "simulation_profile.prof")
    try:
        profiler.dump_stats(profile_output_filename)
        print(f"Full profiling stats saved to {profile_output_filename}")
    except Exception as e_prof_save:
        print(f"Error saving full profile stats: {e_prof_save}")

import os, subprocess, glob, shutil, sys
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", '--processes', type = int, help = "number of parallel processes", default = 30)
    parser.add_argument("-n", '--num_iter', type = int, help = 'number of iterations before restart, needed to tackle memory leakage problems', default = 10000)
    parser.add_argument("-o", '--outfolder', type = str, help = "Path to outfolder where pi_g_results are being stored.")
    parser.add_argument("-m", '--max_files_per_folder', type = int, help = "Maximum number of files allowed to be written to a folder", default = 500_000)

    args = parser.parse_args()

    procs = args.processes
    iters = args.num_iter
    outfolder = args.outfolder
    max_files_per_folder = args.max_files_per_folder

    # Initiate counters
    files_computed  = 0
    first_index     = 0
    folder_idx = 1

    with open("data/total_comparisons.txt", 'r') as fh:
        total_comparisons = int(fh.read().splitlines()[0])

    while files_computed < total_comparisons:
        last_index = first_index + iters

        # This command shall accept a starting index and number of files to be processed in a run.
        # The starting index fot the first time around is zero, next time it is the stop index of the last run. 

        # Working on it in another tab
        cmd = f'python3 scripts/pi_g_MCMC_inference.py -d ./data -p {procs} -n {iters} -o {outfolder} -s {first_index} -e {last_index}'
        
        try:
            subprocess.call(cmd, shell=True)
        except BrokenPipeError:
            pass
        except RuntimeError:
            pass
        
        # We now find the index of the final file
        files_computed = last_index
        first_index    = last_index
        
        files_in_folder = len(glob.glob(os.path.join(outfolder, 'pi_g_results/*.txt')))
        
        # Move results folder if maximum has been reached
        if files_in_folder > max_files_per_folder:
            new_foldername = f"completed_pigs_{folder_idx}"
            shutil.move(os.path.join(outfolder, 'pi_g_results'), os.path.join(outfolder, new_foldername))
            os.mkdir(os.path.join(outfolder, 'pi_g_results'))
            
            folder_idx += 1


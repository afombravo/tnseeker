import multiprocessing
from datetime import datetime
from colorama import Fore
import csv

def cpu():
    c = multiprocessing.cpu_count()
    if c >= 2:
        c -= 1
    return c

def colourful_errors(warning_type,error):
    warning_colour = Fore.GREEN
    if warning_type == "FATAL":
        warning_colour = Fore.RED
    elif warning_type == "WARNING":
        warning_colour = Fore.YELLOW

    print(f"{Fore.BLUE} {datetime.now().strftime('%c')}{Fore.RESET} [{warning_colour}{warning_type}{Fore.RESET}] {error}")

def csv_writer(output_file_path,output_file):
    with open(output_file_path, "w", newline='') as output:
        writer = csv.writer(output)
        writer.writerows(output_file)
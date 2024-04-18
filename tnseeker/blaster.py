import subprocess

def subprocess_cmd(command):
    try:
        return subprocess.check_output(command)
    except subprocess.CalledProcessError as e:
        return e.output.decode()

fasta = "/mnt/c/Users/Afonso/Downloads/Carlos/BU_ATCC8492_annotations/blast/BU_ATCC8492_annotations.fasta"
true_positives = "/mnt/c/Users/Afonso/Downloads/Carlos/BU_ATCC8492_annotations/blast/true_positives.fasta"
true_negatives = "/mnt/c/Users/Afonso/Downloads/Carlos/BU_ATCC8492_annotations/blast/true_negatives.fasta"
blast_db = "/mnt/c/Users/Afonso/Downloads/Carlos/BU_ATCC8492_annotations/blast"
out_pos = f"{blast_db}/blast_pos.txt"
out_neg = f"{blast_db}/blast_neg.txt"

send = ["makeblastdb",
        "-in",
        fasta,
        "-dbtype",
        "nucl"]

subprocess_cmd(send)

for file in [[true_positives,out_pos],
              [true_negatives,out_neg]]:
    
    send = ["tblastn",
            "-query",
            file[0],
            "-db",
            fasta,
            "-evalue",
            "0.00001",
            "-out",
            file[1],
            "-outfmt",
            "6"]
    
    subprocess_cmd(send)

    found_true = set()
    with open(file[1]) as current:
        for line in current:
            line = line.split("\t")
            found_true.add(line[0])
            
    print(len(found_true))
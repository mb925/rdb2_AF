import subprocess
from Bio import SeqIO
from io import StringIO

class Blast:

    def __init__(self, exe_path: str = 'blastdbcmd'):
        # Store path to executable
        self.exe_path = exe_path

    def __call__(self, database_path: str, uniprot_entry: str, cwd_path: str = '.'):
        # Define executable arguments
        exe_args = [self.exe_path]
        exe_args = [*exe_args, '-db', database_path]
        exe_args = [*exe_args, '-entry', uniprot_entry]
        # Debug
        print('Issued Blast command is ' + ' '.join(exe_args))
        # Run process, retrieve output
        # NOTE Encoding is necessary to cast result to string
        return subprocess\
            .run(exe_args, cwd=cwd_path, capture_output=True, encoding='utf-8')\
            .stdout

    @staticmethod
    def get_sequence(text: str) -> tuple:
        fasta = StringIO(text)
        record = SeqIO.parse(fasta, "fasta")
        title = ''
        seq = ''
        for rec in record:
            title = rec.description
            seq = str(rec.seq)
        fasta.close()

        return title, seq

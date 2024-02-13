import os
import pandas as pd
import src.exceptions as e

class SampleInfo(object):
    """Object that takes a minimal input from user
    and fins the paths to tumour, UMI and normal bam files.

    Requires a samplesheet with columns "sample" and "normal"
    and a directory with bam files for tumour, UMI and normal samples.

    """

    def __init__(self, samplesheet, tumour_bam_dir: str, umi_bam_dir: str, normal_bam_dir: str):
        super(SampleInfo, self).__init__()
        self.samplesheet = samplesheet
        self.tumour_bam_dir = tumour_bam_dir
        self.umi_bam_dir = umi_bam_dir
        self.normal_bam_dir = normal_bam_dir
        self.samples = samplesheet["sample"].tolist()
        # find paths to bam files
        self.samplepaths = self.find_bam_files(self.samplesheet, self.tumour_bam_dir)
        self.sample_uncons = self.find_bam_files(self.samplesheet, self.umi_bam_dir)
        self.normals = self.find_bam_files(self.samplesheet, self.normal_bam_dir, sample_col="normal")
        # sampleID : NormalID dict
        self.normal_ids = self.create_normaldicts(self.samplesheet)

    def absoluteFilePaths(self, directory):
        """Yields all file paths in a directory,
        as a generator object."""
        for dirpath,_,filenames in os.walk(directory):
            for f in filenames:
                yield os.path.abspath(os.path.join(dirpath, f))

    def find_bam_files(self, samplesheet: pd.DataFrame, targ_dir: str, sample_col: str = "sample") -> dict:
        """Find all BAM files in a directory for samples found in
        samplesheet. Outputs a dict of SampleID:Path pairs.
        
        """
        # create a list that has all of the file paths in target directories
        bams_all = [path for path in self.absoluteFilePaths(targ_dir) if path.endswith(".bam")]

        out_dict = {}

        for path in bams_all:
            for i, row in samplesheet.iterrows():
                if row[sample_col] in path:
                    out_dict[row['sample']] = path


        # check that the number of keys in out_dict equals the number of samples in the samplesheet
        if len(out_dict.keys()) != len(samplesheet[sample_col].unique().tolist()):
            raise e.FunkyNumberofBams(sample_col, len(out_dict.keys()))

        return out_dict

    def create_normaldicts(self, samplesheet: pd.DataFrame) -> dict:

        normal_ids = {}

        for i, row in samplesheet.iterrows():
            normal_ids[row["sample"]] = row['normal']

        return normal_ids
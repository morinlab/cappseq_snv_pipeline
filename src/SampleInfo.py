import os
import pandas as pd
import warnings
warnings.simplefilter("always")

class SampleInfo(object):
    """Object that takes a minimal input from user
    and fins the paths to tumour, UMI and normal bam files.

    Requires a samplesheet with columns "sample" and "normal"
    and a directory with bam files for tumour, UMI and normal samples.

    Args:
        samplesheet (pd.DataFrame): A pandas dataframe with columns "sample" and "normal".
        tumour_bam_dir (str): The directory containing the tumour bam files.
        umi_bam_dir (str): The directory containing the UMI bam files.
        normal_bam_dir (str): The directory containing the normal bam files.

    Attributes:
    
        samplesheet (pd.DataFrame): A pandas dataframe with columns "sample" and "normal".
        tumour_bam_dir (str): The directory containing the tumour bam files.
        umi_bam_dir (str): The directory containing the UMI bam files.
        normal_bam_dir (str): The directory containing the normal bam files.
        samples (list): A list of sample IDs.
        samplepaths (dict): A dictionary with sample IDs as keys and tumour bam file paths as values.
        sample_uncons (dict): A dictionary with sample IDs as keys and UMI bam file paths as values.
        normals (dict): A dictionary with sample IDs as keys and normal bam file paths as values.
        normal_ids (dict): A dictionary with sample IDs as keys and normal IDs as values.

    """

    def __init__(self, samplesheet, tumour_bam_dir: str, umi_bam_dir: str, normal_bam_dir: str):
        super(SampleInfo, self).__init__()
        self.samplesheet = samplesheet
        self.tumour_bam_dir = tumour_bam_dir
        self.umi_bam_dir = umi_bam_dir
        self.normal_bam_dir = normal_bam_dir

        # make empty log file to be filled with missing samples
        self.make_empty_log("missing_bams_log.txt")
        # find paths to bam files
        self.samplepaths = self.find_bam_files(self.samplesheet, self.tumour_bam_dir)
        self.sample_uncons = self.find_bam_files(self.samplesheet, self.umi_bam_dir)
        self.normals = self.find_bam_files(self.samplesheet, self.normal_bam_dir, sample_col="normal")
        # make sample list from samplepaths keys
        self.samples = list(self.samplepaths.keys())
        # sampleID : NormalID dict
        self.normal_ids = self.create_normaldicts(self.samplesheet)
        self.print_warnings()

    def absoluteFilePaths(self, directory):
        """Yields all file paths in a directory,
        as a generator object
        
        Args:
            directory (str): The directory to search for files."""
        for dirpath,_,filenames in os.walk(directory):
            for f in filenames:
                yield os.path.abspath(os.path.join(dirpath, f))
    
    def make_empty_log(self, log_name: str) -> None:
        """Create an empty log file.
        
        Args:
            log_name (str): The name of the log file to create."""
        with open(log_name, "w") as log:
            log.write("")

    def find_bam_files(self, samplesheet: pd.DataFrame, targ_dir: str, sample_col: str = "sample") -> dict:
        """Find all BAM files in a directory for samples found in
        samplesheet. Outputs a dict of SampleID:Path pairs.

        Also appends missing samples to a log file and warns the user about them in the console.

        Args:
            samplesheet (pd.DataFrame): A pandas dataframe with a column for sample IDs.
            targ_dir (str): The directory to search for BAM files.
            sample_col (str): The name of the column in samplesheet that contains the sample IDs.
        
        Returns:
            dict: A dictionary with sample IDs as keys and BAM file paths as values.
        
        """
        # create a list that has all of the file paths in target directories
        bams_all = [path for path in self.absoluteFilePaths(targ_dir) if path.endswith(".bam")]

        out_dict = {}
        missing_samples = []
        # start a log file to log missing samples, make it easier for user
        with open("missing_bams_log.txt", "a") as log:
            for i, row in samplesheet.iterrows():
                for path in bams_all:
                    if row[sample_col] in path:
                        out_dict[row['sample']] = path
                # test if a dict key exists for each sample in the samplesheet
                
                if row['sample'] not in out_dict.keys():
                    missing_samples.append(row[sample_col])
                    mesg = f"{sample_col} called {row[sample_col]} not found in {targ_dir}"
                    log.write(mesg + "\n")
                    warnings.warn(mesg + '\033[93m')

            # check that the number of unique values in out_dict equals the number of samples in the samplesheet
            # aka one path per sample and normal
            values = set(out_dict.values())
            values = list(values)
            if len(values) != len(samplesheet[sample_col].unique().tolist()):
                mesg2 = f"""For the samples in the {sample_col} column I found {len(out_dict.keys())} bams, but
                    I was expecting {len(samplesheet[sample_col].unique().tolist())}.""" + "\N{loudly crying face}"
                log.write(mesg2 + "\n")
                warnings.warn(mesg2)

            
        
        # write missing samples to tsv
        if len(missing_samples) > 0 and sample_col == "sample":
            pd.DataFrame(missing_samples, columns=["missing_samples"]).to_csv("missing_samples.tsv", sep="\t", index=False)
        elif len(missing_samples) > 0 and sample_col == "normal":
            pd.DataFrame(missing_samples, columns=["missing_normals"]).to_csv("missing_normals.tsv", sep="\t", index=False)

        return out_dict
    
    def print_warnings(self) -> None:
        print("""All warnings about missing samples were written to a log file called 'missing_bams_log.txt'.
        Lists of missing sample and normal IDs were written to 'missing_samples.tsv' and
              'missing_normals.tsv'.""")

    def create_normaldicts(self, samplesheet: pd.DataFrame) -> dict:
        """Create a dictionary with sample IDs as keys and normal IDs as values.

        Args:   
            samplesheet (pd.DataFrame): A pandas dataframe with columns "sample" and "normal".

        Returns:
            dict: A dictionary with sample IDs as keys and normal IDs as values.
        """

        normal_ids = {}

        for i, row in samplesheet.iterrows():
            normal_ids[row["sample"]] = row['normal']

        return normal_ids
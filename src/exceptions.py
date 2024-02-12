# make custom exception called FunkyNumberofBAms
class FunkyNumberofBams(Exception):
    def __init__(self,sample, bam_paths):
        super().__init__()
        self.sample = sample,
        self.bam_paths = bam_paths
    def __str__(self) -> str:
        return (f"Sample {self.sample} had {self.bam_paths} bams, I was expecting 1 only.")
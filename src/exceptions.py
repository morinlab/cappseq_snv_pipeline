# make custom exception called FunkyNumberofBAms
class FunkyNumberofBams(Exception):
    def __init__(self,sample, bam_paths, sample_count):
        super().__init__()
        self.sample = sample
        self.bam_paths = bam_paths
        self.sample_count = sample_count
    def __str__(self) -> str:
        return (f"""For the samples in the {self.sample} column I found {self.bam_paths} bams, but
                I was expecting {self.sample_count}.""" + "\N{loudly crying face}")
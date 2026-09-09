WDL script that converts a VCF into a MAPLE-formatted diff file, for use with [UShER](https://github.com/yatisht/usher) (specifically `usher sampled diff`). If you are calling *Mycobacterium* samples with [clockwork](https://github.com/iqbal-lab-org/clockwork), and want to put them on a tree with UShER, you are in the right place. If you are working with viruses such as SARS-CoV-2 you probably are using a different usher call and don't need this.

This is mostly just a wrapper of the VCF->diff conversion Python script written by Lily Karim: https://github.com/lilymaryam/parsevcf

import sib_regrid
import sib_monthly_combiner
import datetime

if __name__ == "__main__":
    print datetime.datetime.now()
    sib_regrid.LRU_paper_main(year=2008, month=7)
    sib_regrid.LRU_paper_main(year=2008, month=8)
    sib_monthly_combiner.LRU_paper_main()
    print datetime.datetime.now()

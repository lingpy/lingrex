all:
	@python general.py
	@python predict.py data/burmish-240-8.tsv -r 0.75 --runs 2 
	@python predict.py data/chinese-623-14.tsv -r 0.75 --runs 2
	@python predict.py data/polynesian-210-10.tsv -r 0.75 --runs 2
	@python predict.py data/japanese-200-10.tsv -c crossid -r 0.75 --runs 2




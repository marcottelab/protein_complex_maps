

cat elutions.txt | parallel --memfree (size) 'python /home/kdrew/scripts/blake_complexes/score.py --alkj {} --lskdj {.}_poissoncorr.txt --optlaskdj poisson'

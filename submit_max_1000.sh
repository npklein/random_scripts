for dr in $(seq 0 0.1 1);
do
  echo DR${dr}
  cd DR${dr}
  for i in {1..22};
  do
    echo chr${i}
    cd chr${i}
    for f in *sh;
    do
      while :
      do
        numberOfJobs=$(squeue -u umcg-ndeklein | wc -l)
        if [ "$numberOfJobs" -lt "1002" ];
        then
           break
        fi
        echo "$numberOfJobs of my jobs in queue, waiting to submit more. Sleeping..."
        sleep 10
      done
      if [ ! -f ${f%sh}err ];
      then
        sbatch $f
      else
        echo "$f already has .err file, skip"
      fi
   done
   cd ../
  done
  cd ../
done

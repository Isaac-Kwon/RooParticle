
echo $1 $2 $3
echo RUN $1 $2 : $3 times

njob=10

for (( i = 0; i < $3; )); do
  for (( j = 0; j < $njob-1; i++&&j++)); do
    $1 $2 $i &
    sleep 2
  done
  if [ ${i} -eq $3 ];then
    break
  fi
  $1 $2 $i
  i=$((i+=1))
  # echo outside
done

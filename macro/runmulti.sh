
echo $1 $2 $3

echo RUN $1 $2 times

for (( i = 0; i < $3; i++ )); do
  echo $1 $2 $i
done

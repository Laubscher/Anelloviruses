#! /bin/bash
#sequence in 1 line and header cut after "."
>1l_$1
for ligne in $(cat $1 | cut -d"." -f 1)
do
if [ "${ligne::1}" = ">" ]
then
echo ${ligne} >> 1l_$1
else
echo -n ${ligne} >> 1l_$1
fi
done
sed -i s"/>/\n>/" 1l_$1

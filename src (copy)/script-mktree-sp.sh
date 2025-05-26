for j in $(ls); do cat $j/*.fasta | sed -E s/0$/0$"-$j-"/g | sed -E s/1$/1$"-$j-"/g | sed -E s/2$/2$"-$j-"/g | sed -E s/3$/3$"-$j-"/g | sed -E s/4$/4$"-$j-"/g | sed -E s/5$/5$"-$j-"/g | sed -E s/6$/6$"-$j-"/g | sed -E s/7$/7$"-$j-"/g | sed -E s/8$/8$"-$j-"/g | sed -E s/9$/9$"-$j-"/g | sed -E s/-$//g; done



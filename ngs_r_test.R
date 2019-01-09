normal_sample <- runif(10000)
ds_matrix = matrix(c(1,2,3,4,5,6), nrow = 3, ncol=2)

hist(ds_matrix, col= "yellow", xlab='seed position')

vec_a = c(1,2,3,4,5)
names(vec_a) = c('id1', 'id2', 'id3', 'id4', 'id5')
vec_a

rownames(ds_matrix) = c('row1', 'row2', 'row3')
colnames(ds_matrix) = c('col1', 'col2')

vec_name = c('Michael', 'Kevin', 'Steve')
vec_height = c(174, 160, 182)
vec_weight = c(70, 66, 80)
df_student = data.frame(name=vec_name, height=vec_height, weight = vec_weight)
df_student$bmi= df_student$weight / (df_student$height * df_student$height)
df_student$bmi = df_student$bmi* 100^2
df_student <- df_student[, 1:3]
tmp = data.frame(name='sandy', height=163, weight=59)
df_student = rbind(tmp, df_student)
df_student

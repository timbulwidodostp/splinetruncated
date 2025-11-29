# Olah Data Semarang
# WhatsApp : +6285227746673
# IG : @olahdatasemarang_
# Nonparametric Spline Truncated Regression Use splinetruncated With (In) R Software
install.packages("pracma")
install.packages("car")
splinetruncated = read.csv("https://raw.githubusercontent.com/timbulwidodostp/splinetruncated/main/splinetruncated/splinetruncated.csv",sep = ";")
# Estimation Nonparametric Spline Truncated Regression Use splinetruncated With (In) R Software
source('https://raw.githubusercontent.com/timbulwidodostp/splinetruncated/main/splinetruncated/splinetruncated.R')
x = splinetruncated[, -c(1,2)]
y = splinetruncated[, 2]
model = spline.truncated.linier.multivariabel(x, y, b = 30, taraf.alpha = 0.05)
# Nonparametric Spline Truncated Regression Use splinetruncated With (In) R Software
# Olah Data Semarang
# WhatsApp : +6285227746673
# IG : @olahdatasemarang_
# Finished
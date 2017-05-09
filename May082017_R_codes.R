3 + 5
12/7
weight_kg = 55
weight_kg

weight_kg = 57.5

double = 2.2
double * weight_kg

round(126.5)
sqrt(126.5)

round(11.24722, digits=2)
round(sqrt(126.5), digits =2)
round(weight_kg)

weight_g<-c(50, 60, 65, 82)
weight_g

animals <-c("mouse", "rat", "dog")
animals
length(animals)

class(animals)
weight_g
class(weight_g)
str(weight_g)
str(animals)
animals<-c(animals, c("cat", "dog"))
animals
str(animals)

animals
animals[2]
animals[c(1, 5)]
animals[c(1:3)]

animals
animals[c(TRUE, TRUE, FALSE, FALSE, TRUE)]

weight_g
weight_g<-c(21, 34, 39, 54, 55)
weight_g
weight_g[c(TRUE, FALSE, TRUE, TRUE, FALSE)]

weight_g > 50
weight_g[weight_g >50]
weight_g[weight_g >50 | weight_g <30]
weight_g
weight_g[weight_g>=30 & weight_g == 21]
weight_g[weight_g>=30 | weight_g == 21]

animals[animals == "cat" | animals == "rat"]
animals[animals %in% c("mouse", "cat")]
animals[animals %in% c("mouse", "cat", "cammel")]

small_list <- c("mouse", "cat", "cammel")
animals[animals %in% small_list]

heights <- c(2, 4, 4, NA, 6)
heights
mean(heights)
max(heights)
mean(heights, na.rm=TRUE)
max(heights, na.rm=T)

heights
heights[is.na(heights)]
heights[! is.na(heights)]

na.omit(heights)

lengths<-c(10, 2, 4, NA, 18, NA, 20)
median(lengths)
median(lengths[! is.na(lengths)])


download.file("https://ndownloader.figshare.com/files/2292169","portal_data_joined.csv")
surveys <- read.csv("portal_data_joined.csv")
head(surveys)
tail(surveys)
dim(surveys)
class(surveys)
str(surveys)

nrow(surveys)
ncol(surveys)
names(surveys)

head(surveys)
rownames(surveys)
rownames(surveys)[1:5]
names(surveys)[1:5]

summary(surveys)
names(surveys)
head(surveys$species_id)
str(surveys)

length(unique(surveys$species_id))

surveys[1]
dim(surveys)
surveys[1, 1]
surveys[2, 1]
surveys[, 1]
surveys[1, ]
length(surveys[1:3, 7])
head_surveys <-surveys[1:6, ]
dim(head_surveys)
new_surveys<-surveys[, -1]
dim(new_surveys)
dim(surveys)
new_surveys<-surveys[, -5]
new_surveys<-surveys[, -c(1:5)]
new_surveys<-surveys[, -c(1:5, 10)]

names(surveys)
head(surveys[, "species_id"])
head(surveys$species_id)
surveys["species_id"]
head(surveys[["species_id"]])

str(surveys)

sex <- factor(c("male", "female", "female", "male"))
sex
levels(sex)
nlevels(sex)
sex
sex<-factor(sex, levels=c("male", "female"))
sex

as.character(sex)
str(as.character(sex))

f<-factor(c(1990, 1983, 1977, 1998, 1990))
f
as.numeric(f)
as.numeric(as.character(f))

sex
plot(surveys$sex)

sex<-surveys$sex
head(sex)
levels(sex)
levels(sex)[1]<-"missing"
levels(sex)
head(sex)

levels(sex)
levels(sex)<-c("missing", "female", "male")
head(sex)
sex<-factor(sex, levels=c("female", "male", "missing"))
plot(sex)

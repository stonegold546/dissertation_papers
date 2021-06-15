# Load and describe dataset ----
dat <- read.csv(paste0(data.location, "mpi.csv"))
colnames(dat)

# order by country
dat <- dat[order(dat$Country), ]

# Columns 7-12 are living standards data
describe(dat)

# Bring in CO2 data ----
# co2.dat <- read.csv(paste0(data.location, "co2.csv"))
# colnames(co2.dat)
# colnames(co2.dat)[1] <- "Country"
# co2.dat$ln_co2 <- log(co2.dat$YR2014)

# Merge data to CO2 emissions ----
# dat <- merge(dat, co2.dat[, c("Country", "ln_co2")], all.x = TRUE)
# rm(co2.dat)

# Make data long format ----
X <- reshape(
  dat, direction = "long", varying = 7:12, v.names = "response", timevar = "item")
# Add item names
X$item_n <- factor(X$item, levels = 1:6, labels = gsub(".", " ", colnames(dat)[7:12], fixed = TRUE))
X$item_t <- paste0("x", X$item)  # Add item id: x1:x6
X <- X[order(X$id), ]  # Sort by country
X$response <- X$response / 100  # Rescale percentages to 0-1 interval
X <- X[!is.na(X$response), ]  # Drop missing data, 3 cases
X$response.01 <- X$response  # Duplicate responses
# Check minimum and maximum values
aggregate(response.01 ~ item_n, X, function (x) c(min(x), max(x)))
# Rescale cooking fuel and electricity using transformation from
# Smithson & Verkuilen, 2006
X$response.01[X$item_n == "Cooking fuel"] <-
  (X$response.01[X$item_n == "Cooking fuel"] * (sum(X$item_n == "Cooking fuel") - 1) + .5) /
  sum(X$item_n == "Cooking fuel")
X$response.01[X$item_n == "Electricity"] <-
  (X$response.01[X$item_n == "Electricity"] * (sum(X$item_n == "Electricity") - 1) + .5) /
  sum(X$item_n == "Electricity")
# Check minimum and maximum values again
aggregate(response.01 ~ item_n, X, function (x) c(min(x), max(x)))

# Transform to probits ----
X$response.l <- -qnorm(X$response.01)

# Extend original dataset with 0-1 interval data ----
Xw <- reshape(X[, c("Country", "item", "response.l")], direction = "wide",
              timevar = "item", v.names = "response.l", idvar = "Country")
Xw
dat <- merge(dat, Xw)
colnames(dat)[15:20] <- paste0("x", 1:6)
rm(Xw)


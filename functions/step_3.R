library(ggplot2)
# For Page 7
# 3.1
## Prevalence for those aged between 40 and 44
obs_data = list(pos=c(55, 145, 183, 164, 212, 299, 222, 190, 133),
                tot=c(307, 797, 829, 633, 598, 642, 481, 439, 366))
obs_data$pos[6] / obs_data$tot[6] # prevalence: 0.4657321
# CI: 0.446 - 0.485
obs_data$pos[6] / obs_data$tot[6] + c(-1,1) * sqrt(obs_data$pos[6] / obs_data$tot[6] * (1 - obs_data$pos[6] / obs_data$tot[6]) / obs_data$tot[6])

## Prevalence for those aged under 25
sum(obs_data$pos[1:2]) / sum(obs_data$tot[1:2]) # prevalence: 0.1811594
# CI: 0.170 - 0.193
sum(obs_data$pos[1:2]) / sum(obs_data$tot[1:2]) + c(-1,1) * sqrt(sum(obs_data$pos[1:2]) / sum(obs_data$tot[1:2]) * (1 - sum(obs_data$pos[1:2]) / sum(obs_data$tot[1:2])) / sum(obs_data$tot[1:2]))

## For mean years of survival, I don't know
# Recidivism Probability
rho <- c(
  0.5472816,
  0.6104789,
  0.5816654,
  0.5151572,
  0.4358284,
  0.4178032,
  0.4031776,
  0.4029714,
  0.4315638
)
rearrest = cbind(t=seq(0,10,0.1),as.data.frame(do.call(cbind,lapply(rho,function(x) exp(-x*seq(0,10,0.1))))))
colnames(rearrest) = c('t','15-19','20-24','25-29','30-34','35-39','40-44','45-49','50-54','55+')
# Draw all 11 lines in one plot with scale bar
rearrest %>% 
  ggplot() + 
  geom_line(aes(x = t, y = `15-19`, color = '15-19')) + 
  geom_line(aes(x = t, y = `20-24`, color = '20-24')) + 
  geom_line(aes(x = t, y = `25-29`, color = '25-29')) + 
  geom_line(aes(x = t, y = `30-34`, color = '30-34')) + 
  geom_line(aes(x = t, y = `35-39`, color = '35-39')) + 
  geom_line(aes(x = t, y = `40-44`, color = '40-44')) + 
  geom_line(aes(x = t, y = `45-49`, color = '45-49')) + 
  geom_line(aes(x = t, y = `50-54`, color = '50-54')) + 
  geom_line(aes(x = t, y = `55+`, color = '55+')) + 
  theme_minimal() + 
  labs(title = 'Recidivism Probability by Age Group', x = 'Years', y = 'Recidivism Probability') +
  scale_color_manual(values = c('15-19' = 'red', '20-24' = 'blue', '25-29' = 'green', 
                                '30-34' = 'purple', '35-39' = 'orange', '40-44' = 'black', 
                                '45-49' = 'brown', '50-54' = 'pink', '55+' = 'grey')) + 
  theme(legend.position = 'right',
        plot.title = element_text(hjust = 0.5))
# Median sojourn time for 35+
median(1 / rho[5:9]) #2.393471

# 3.2
revise = matrix(0, nrow = 100,ncol = 48)
for (iter in 1:100) {
  # Load each iteration's result sampled from posterior distribution
  TAB = readRDS(paste0('output/BI_normal_SVR/sim_data_normal_SVR_',iter,'.rds'))
  
  TAB$t = TAB$t/360 # convert to 'years'
  
  strategies = read.csv("settings/Final_strategies.csv")
  strategies$pD = strategies$pD/12;   # Original is *loosely the annual proportion
  strategies$pJ1 = strategies$pJ1/12; # "
  strategies$pJ2 = strategies$pJ2/12; # "
  strategies$pJ3 = strategies$pJ3/12; # "
  strategies$pX = strategies$pX/12;   # "
  
  # Extract needed results
  D <- cbind(TAB[,1:2], TAB[ , grepl( "D" , names( TAB ) ) ])
  J <- cbind(TAB[,1:2], TAB[ , grepl( "J" , names( TAB ) ) ])
  X <- cbind(TAB[,1:2], TAB[ , grepl( "X" , names( TAB ) ) ])
  
  C4 <- cbind(TAB[,1:2], TAB[ , grepl( "C4" , names( TAB ) ) ])
  C5 <- cbind(TAB[,1:2], TAB[ , grepl( "C5" , names( TAB ) ) ])
  C6 <- cbind(TAB[,1:2], TAB[ , grepl( "C6" , names( TAB ) ) ])
  
  ## Targeting the current PWID and treating 15% of those with mild to compensated 
  ## cirrhosis would potentially lower the total number of infections to under 10,000 
  ## by year 50, with counts among PWID or those in prisons drop to 118 or lower.
  revise[iter,1] = tail(rowSums(D[which(TAB$Strategy==19),12:74]) + rowSums(J[which(TAB$Strategy==19),12:74]),1)
  
  ## with slightly more cases to be averted among former PWID in the community (by 15259)
  ## and fewer reductions among current PWID (by 8512) by year 50
  revise[iter,2] = tail(rowSums(X[which(TAB$Strategy==1),12:74])-rowSums(X[which(TAB$Strategy==2),12:74]),1)
  revise[iter,3] = tail(rowSums(D[which(TAB$Strategy==1),12:74])-rowSums(D[which(TAB$Strategy==2),12:74]),1)
  
  ## all incarcerated ex-PWID, and 12% of the formerly incarcerated ex-PWID diagnosed 
  ## with mild to compensated cirrhosis—could reduce the number of HCV infection to 
  ## fewer than 618 within the next 30 years
  revise[iter,4] = cbind(TAB$t[which(TAB$Strategy==37)],rowSums(D[which(TAB$Strategy==37),12:74]+J[which(TAB$Strategy==37),12:74]+X[which(TAB$Strategy==37),12:74]))[1+12*30,2]
  
  # For Page 8
  ## The combined strategy remains the most effective and efficient in averting 
  ## HCV-related complications, substantially lowering the annual number of new 
  ## cases of decompensated cirrhosis and hepatocellular carcinoma to fewer than 
  ## 2697 and 402, respectively, within the first 20 years of intervention
  revise[iter,5] = cbind(TAB$t[which(TAB$Strategy==38)],rowSums(C4[which(TAB$Strategy==38),3:11]))[1+12*20,2]
  revise[iter,6] = cbind(TAB$t[which(TAB$Strategy==38)],rowSums(C5[which(TAB$Strategy==38),3:11]))[1+12*20,2]
  
  # # of strategies
  str = 1
  
  # For Page 9
  # 3.3
  # Table 3
  # Add strategy 35 to check what it looks like when proportion of Ex-PWID in community is about 10% in combined strategy
  for (row in c(1,19,2,34,38,35)) {
    df = TAB %>% filter(Strategy == row, t == 15)
    
    # Prevalence
    # All
    revise[iter,7 * str] = 100 - 100*rowSums(select(df,grep("u", names(df), value = TRUE))) / rowSums(select(df, matches('^[DJX]')))
    
    # PWID
    revise[iter,1 + 7 * str] = 100 - 100*rowSums(select(df,grep("^Du[0-9]", names(df), value = TRUE))) / rowSums(select(df,grep("D", names(df), value = TRUE)))
    
    # In Prison
    revise[iter,2 + 7 * str] = 100 - 100*rowSums(select(df,grep("^Ju[0-9]", names(df), value = TRUE))) / rowSums(select(df,grep("J", names(df), value = TRUE)))
    
    # Ex-PWID
    revise[iter,3 + 7 * str] = 100 - 100*rowSums(select(df,grep("^Xu[0-9]", names(df), value = TRUE))) / rowSums(select(df,grep("X", names(df), value = TRUE)))
    
    # DC number
    revise[iter,4 + 7 * str] = df %>% select(grep("^C4[0-9]", names(df), value = TRUE)) %>% rowSums()
    
    # HC number
    revise[iter,5 + 7 * str] = df %>% select(grep("^C5[0-9]", names(df), value = TRUE)) %>% rowSums()
    
    # Treated number
    revise[iter,6 + 7 * str] = df$ctx
    
    str = str + 1
  }
  
  cat('Iteration',iter,'\n')
}

saveRDS(revise,'output/revise_matrix.rds')

# Appendix
# Figure S2
# Load the grid package
library(grid)

# Load posterior samples
pred_prev <- readRDS('output/mcmc/mcmc_bci.rds')
pred_prev_median <- apply(pred_prev, 2, median)
pred_prev_lower <- apply(pred_prev, 2, quantile, probs = 0.025)
pred_prev_upper <- apply(pred_prev, 2, quantile, probs = 0.975)

# Hard-coded observations
obs_data <- list(
  pos = c(55, 145, 183, 164, 212, 299, 222, 190, 133),
  tot = c(307, 797, 829, 633, 598, 642, 481, 439, 366)
)
age_grp <- c("15-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-49", "50-54", "55+")
percentages <- 100 * c(
  0.1791531, 0.1819322, 0.2207479, 0.2590837, 0.3545151, 0.4657321,
  0.4615385, 0.4328018, 0.3633880
)
errors <- 100 * sqrt(obs_data$pos / obs_data$tot * (1 - obs_data$pos / obs_data$tot) / obs_data$tot)

# Fitted prevalence and errors from model
fitted_percentages <- 100 * pred_prev_median
fitted_upper_errors <- 100 * (pred_prev_upper - pred_prev_median)
fitted_lower_errors <- 100 * (pred_prev_median - pred_prev_lower)

# Original colors
colors <- c(
  "#FFCB47", # Light Yellow-Orange
  "#FFB65C", # Soft Orange
  "#FF9A1F", # Vivid Orange
  "#FB804B", # Coral Orange
  "#FB6423", # Deep Orange
  "#FF4747", # Bright Red
  "#FF1F1F", # Vivid Red
  "#D0253C", # Dark Red
  "#B8282B" # Deep Brick Red
)
grey_colors <- rev(c("#D3D3D3", "#C0C0C0", "#B0B0B0", "#A9A9A9", "#A0A0A0", "#989898", "#909090", "#888888", "#808080"))

# Fitted colors
fitted_colors <- c(
  "#00D4FF", # Bright Cyan
  "#01BCEA", # Sky Blue
  "#02A5D5", # Medium Cyan-Blue
  "#038DC0", # Deep Aqua Blue
  "#0575AB", # Steel Blue
  "#065D95", # Dark Blue
  "#074680", # Navy Blue
  "#082E6B", # Indigo Blue
  "#091656" # Deep Navy
)

# Create a viewport for the entire plot area
grid.newpage()
pushViewport(viewport(width = 1.0, height = 1.0))

# Y-axis label
grid.text("HCV Prevalence (%)", x = unit(0.03, "npc"), y = unit(0.5, "npc"), rot = 90, gp = gpar(fontsize = 21))

# X-axis label
grid.text("Age Group", x = unit(0.5, "npc"), y = unit(0.1, "npc"), gp = gpar(fontsize = 21))

# Draw the bars with error bars
bar_width <- 0.07 # Width of bars
error_bar_width <- bar_width / 3 # Width of error bars
scale <- 1.8 # Scale factor for y-axis

# Y-axis ticks and labels
for (i in seq(10, 50, by = 10)) {
  grid.text(i, x = unit(0.09, "npc"), y = unit(0.2 - 0.1 * scale + i / 100 * scale, "npc"), just = "right", gp = gpar(fontsize = 19))
  grid.lines(x = unit(c(0.1, 0.11), "npc"), y = unit(0.2 - 0.1 * scale + i / 100 * scale, "npc"))
}

# Add axis line
# X-axis
grid.lines(x = unit(c(0.11, 0.97), "npc"), y = unit(0.2, "npc"))
# Y-axis
grid.lines(x = unit(0.11, "npc"), y = unit(c(0.2, 0.95), "npc"))

for (i in 1:length(age_grp)) {
  x_center <- 0.111 - 2 * error_bar_width + i * 4 * error_bar_width

  # Coordinates
  y_obs <- 0.2 - 0.1 * scale + scale * (percentages[i]) / 100
  y_fit <- 0.2 - 0.1 * scale + scale * (fitted_percentages[i]) / 100
  y_obs_lo <- 0.2 - 0.1 * scale + scale * (percentages[i] - errors[i]) / 100
  y_obs_hi <- 0.2 - 0.1 * scale + scale * (percentages[i] + errors[i]) / 100
  y_fit_lo <- 0.2 - 0.1 * scale + scale * (fitted_percentages[i] - fitted_lower_errors[i]) / 100
  y_fit_hi <- 0.2 - 0.1 * scale + scale * (fitted_percentages[i] + fitted_upper_errors[i]) / 100

  # ----- Draw vertical error bars -----
  grid.segments(
    x0 = unit(x_center - 0.015, "npc"), x1 = unit(x_center - 0.015, "npc"),
    y0 = unit(y_obs_lo, "npc"), y1 = unit(y_obs_hi, "npc"),
    gp = gpar(col = colors[4], lwd = 5)
  )
  grid.segments(
    x0 = unit(x_center + 0.015, "npc"), x1 = unit(x_center + 0.015, "npc"),
    y0 = unit(y_fit_lo, "npc"), y1 = unit(y_fit_hi, "npc"),
    gp = gpar(col = fitted_colors[4], lwd = 5)
  )

  # ----- Draw dots for means -----
  # Observed
  grid.points(
    x = unit(x_center - 0.015, "npc"),
    y = unit(y_obs, "npc"),
    pch = 21, # solid circle
    size = unit(8, "pt"),
    gp = gpar(fill = colors[4], col = "black")
  )

  # Fitted (shifted slightly to the right)
  grid.points(
    x = unit(x_center + 0.015, "npc"),
    y = unit(y_fit, "npc"),
    pch = 21, # filled circle with border
    size = unit(8, "pt"),
    gp = gpar(fill = fitted_colors[4], col = "black")
  )

  # X-axis tick
  grid.lines(x = unit(x_center + 2 * error_bar_width, "npc"), y = unit(c(0.19, 0.2), "npc"))

  # Label
  grid.text(age_grp[i],
    x = unit(x_center, "npc"),
    y = unit(0.16, "npc"),
    just = "center", gp = gpar(fontsize = 19)
  )
}

# X-axis tick
grid.lines(x = unit(0.11, "npc"), y = unit(c(0.19, 0.2), "npc"))

obs_legend_x <- 0.315
obs_legend_y <- 0.88

fitted_legend_x <- 0.28
fitted_legend_y <- 0.78

grid.text("Observed", x = unit(obs_legend_x + 0.015, "npc"), y = unit(obs_legend_y, "npc"), just = "right", gp = gpar(fontsize = 19))
grid.text("Fitted", x = unit(fitted_legend_x, "npc"), y = unit(fitted_legend_y, "npc"), just = "right", gp = gpar(fontsize = 19))

grid.rect(
  x = unit(0.15, "npc"), y = unit(obs_legend_y, "npc"),
  width = unit(2 * error_bar_width, "npc"), height = unit(error_bar_width, "npc"),
  just = "left",
  gp = gpar(fill = colors[4], col = NA) # no border
)

grid.rect(
  x = unit(0.15, "npc"), y = unit(fitted_legend_y, "npc"),
  width = unit(2 * error_bar_width, "npc"), height = unit(error_bar_width, "npc"),
  just = "left",
  gp = gpar(fill = fitted_colors[4], col = NA) # no border
)


# Update Table 3 and Table S1
revise = readRDS('output/revise_matrix_low_SVR.rds')

for (col in 7:ncol(revise)){
  print(quantile(revise[,col], c(0.025, 0.975)))
  cat('Press Enter to continue','\n')
  readline()
}

TAB = read.csv('settings/tableHIGH_1_wFINALstrategies_yichen_low_SVR.csv')
TAB$t = TAB$t/360 # convert to 'years'

strategies = read.csv("settings/Final_strategies.csv")
strategies$pD = strategies$pD/12;   # Original is *loosely the annual proportion
strategies$pJ1 = strategies$pJ1/12; # "
strategies$pJ2 = strategies$pJ2/12; # "
strategies$pJ3 = strategies$pJ3/12; # "
strategies$pX = strategies$pX/12;   # "
for (row in c(1,19,2,34,38,35)) {
  df = TAB %>% filter(Strategy == row, t == 15)
  
  # Prevalence
  # All
  print(100 - 100*rowSums(select(df,grep("u", names(df), value = TRUE))) / rowSums(select(df, matches('^[DJX]'))))
  
  # PWID
  print(100 - 100*rowSums(select(df,grep("^Du[0-9]", names(df), value = TRUE))) / rowSums(select(df,grep("D", names(df), value = TRUE))))
  
  # In Prison
  print(100 - 100*rowSums(select(df,grep("^Ju[0-9]", names(df), value = TRUE))) / rowSums(select(df,grep("J", names(df), value = TRUE))))
  
  # Ex-PWID
  print(100 - 100*rowSums(select(df,grep("^Xu[0-9]", names(df), value = TRUE))) / rowSums(select(df,grep("X", names(df), value = TRUE))))
  
  # DC number
  print(df %>% select(grep("^C4[0-9]", names(df), value = TRUE)) %>% rowSums())
  
  # HC number
  print(df %>% select(grep("^C5[0-9]", names(df), value = TRUE)) %>% rowSums())
  
  # Treated number
  print(df$ctx)
  
  readline('Press Enter to continue...')
}


# Sensitivity Analysis
for (iter in 1:100) {
  # Load each iteration's result sampled from posterior distribution
  TAB = readRDS(paste0('output/BI_normal_SVR/sim_data_normal_SVR_',iter,'.rds'))
  
  TAB$t = TAB$t/360 # convert to 'years'
  
  strategies = read.csv("settings/Final_strategies.csv")
  # "
  strategies$pD = strategies$pD/12;   # Original is *loosely the annual proportion
  strategies$pJ1 = strategies$pJ1/12; # "
  strategies$pJ2 = strategies$pJ2/12; # "
  strategies$pJ3 = strategies$pJ3/12; # "
  strategies$pX = strategies$pX/12;   # "
}




# Appendix 
## % difference in reduction change of each treatment group
# Mean
df = TAB %>% select(!contains('u')) %>% select(!contains('C'))
df_combined = df %>% rowwise() %>% 
  mutate('.1' = D01+D11+D21+D31+D41+D51+D61+J01+J11+J21+J31+J41+J51+J61+X01+X11+X21+X31+X41+X51+X61,
         '.2' = D02+D12+D22+D32+D42+D52+D62+J02+J12+J22+J32+J42+J52+J62+X02+X12+X22+X32+X42+X52+X62,
         '.3' = D03+D13+D23+D33+D43+D53+D63+J03+J13+J23+J33+J43+J53+J63+X03+X13+X23+X33+X43+X53+X63,
         '.4' = D04+D14+D24+D34+D44+D54+D64+J04+J14+J24+J34+J44+J54+J64+X04+X14+X24+X34+X44+X54+X64,
         '.5' = D05+D15+D25+D35+D45+D55+D65+J05+J15+J25+J35+J45+J55+J65+X05+X15+X25+X35+X45+X55+X65,
         '.6' = D06+D16+D26+D36+D46+D56+D66+J06+J16+J26+J36+J46+J56+J66+X06+X16+X26+X36+X46+X56+X66,
         '.7' = D07+D17+D27+D37+D47+D57+D67+J07+J17+J27+J37+J47+J57+J67+X07+X17+X27+X37+X47+X57+X67,
         '.8' = D08+D18+D28+D38+D48+D58+D68+J08+J18+J28+J38+J48+J58+J68+X08+X18+X28+X38+X48+X58+X68,
         '.9' = D09+D19+D29+D39+D49+D59+D69+J09+J19+J29+J39+J49+J59+J69+X09+X19+X29+X39+X49+X59+X69) %>% 
  select(-c(3:191))

for (row in c(19,2,34,35,38)) {
  temp_df = df_combined %>% filter(t == 50)
  
  # % reduction comparing to baseline group
  print((100 - 100 * (filter(temp_df,Strategy == row) / filter(temp_df,Strategy == 1)))[,-c(1:2)])
  
  cat('Press Enter to continue:\n')
  
  readline()
}


# BCIs
reduction = matrix(0, nrow = 100, ncol = 45)
for (iter in 1:100) {
  # Load each iteration's result sampled from posterior distribution
  TAB = readRDS(paste0('output/BI_normal_SVR/sim_data_normal_SVR_',iter,'.rds'))
  
  TAB$t = TAB$t/360 # convert to 'years'
  
  strategies = read.csv("settings/Final_strategies.csv")
  # "
  strategies$pD = strategies$pD/12;   # Original is *loosely the annual proportion
  strategies$pJ1 = strategies$pJ1/12; # "
  strategies$pJ2 = strategies$pJ2/12; # "
  strategies$pJ3 = strategies$pJ3/12; # "
  strategies$pX = strategies$pX/12;   # "
  
  df = TAB %>% select(!contains('u')) %>% select(!contains('C'))
  df_combined = df %>% rowwise() %>% 
    mutate('.1' = D01+D11+D21+D31+D41+D51+D61+J01+J11+J21+J31+J41+J51+J61+X01+X11+X21+X31+X41+X51+X61,
           '.2' = D02+D12+D22+D32+D42+D52+D62+J02+J12+J22+J32+J42+J52+J62+X02+X12+X22+X32+X42+X52+X62,
           '.3' = D03+D13+D23+D33+D43+D53+D63+J03+J13+J23+J33+J43+J53+J63+X03+X13+X23+X33+X43+X53+X63,
           '.4' = D04+D14+D24+D34+D44+D54+D64+J04+J14+J24+J34+J44+J54+J64+X04+X14+X24+X34+X44+X54+X64,
           '.5' = D05+D15+D25+D35+D45+D55+D65+J05+J15+J25+J35+J45+J55+J65+X05+X15+X25+X35+X45+X55+X65,
           '.6' = D06+D16+D26+D36+D46+D56+D66+J06+J16+J26+J36+J46+J56+J66+X06+X16+X26+X36+X46+X56+X66,
           '.7' = D07+D17+D27+D37+D47+D57+D67+J07+J17+J27+J37+J47+J57+J67+X07+X17+X27+X37+X47+X57+X67,
           '.8' = D08+D18+D28+D38+D48+D58+D68+J08+J18+J28+J38+J48+J58+J68+X08+X18+X28+X38+X48+X58+X68,
           '.9' = D09+D19+D29+D39+D49+D59+D69+J09+J19+J29+J39+J49+J59+J69+X09+X19+X29+X39+X49+X59+X69) %>% 
    select(-c(3:191))
  
  strategy = 0
  
  for (row in c(19,2,34,35,38)) {
    temp_df = df_combined %>% filter(t == 50)
    
    # % reduction comparing to baseline group
    reduction[iter, (9 * strategy + 1:9)] = as.numeric((100 - 100 * (filter(temp_df,Strategy == row) / filter(temp_df,Strategy == 1)))[,-c(1:2)])
    
    strategy = strategy + 1
    #cat('Press Enter to continue:\n')
    
    #readline()
  }
  print(iter)
}
saveRDS(reduction,'output/reduction.rds')


# Section 4
# How many numbers we need to treat
TAB = read.csv('settings/tableHIGH_1_wFINALstrategies_yichen.csv')
TAB$t = TAB$t/360 # convert to 'years'

strategies = read.csv("settings/Final_strategies.csv")
# "
strategies$pD = strategies$pD/12;   # Original is *loosely the annual proportion
strategies$pJ1 = strategies$pJ1/12; # "
strategies$pJ2 = strategies$pJ2/12; # "
strategies$pJ3 = strategies$pJ3/12; # "
strategies$pX = strategies$pX/12;   # "

# I do the average of the first 10 yrs
TAB %>% filter(Strategy==34,t %in% 1:10) %>% pull(ctx) %>% diff %>% mean

TAB %>% filter(Strategy==34,t %in% 11:50) %>% pull(ctx) %>% diff %>% mean






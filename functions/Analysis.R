library(ggplot2)
df_bo = TAB_bo %>% filter(Strategy==1) %>% select(contains('Du')) %>% rowSums()
df_yichen = TAB_yichen %>% filter(Strategy==1) %>% select(contains('Du')) %>% rowSums()

# Plot
df = data.frame(x = 1:length(df_bo),Bo = df_bo, Yichen = df_yichen)
ggplot(df) + 
  geom_line(mapping = aes(x = x,y = Bo, color = 'red')) + 
  geom_line(mapping = aes(x = x,y = Yichen, color = 'blue'))


df_yichen= TAB_yichen %>% 
  filter(Strategy == 2) %>% 
  select(ctx) %>% 
  rowSums() %>%
  round(2)

df_bo = TAB_bo %>% 
  filter(Strategy == 2) %>% 
  select(ctx) %>% 
  rowSums() %>%
  round(2)


## The generation of Table 3 in Manuscript
for (row in c(1,19,2,34,38)) {
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

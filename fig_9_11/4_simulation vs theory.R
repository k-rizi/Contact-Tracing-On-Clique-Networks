
# test_cpp ----------------------------------------------------------------

# Install the Rcpp package if you haven't already
if (!requireNamespace("Rcpp", quietly = TRUE)) {
  install.packages("Rcpp")
}

# Include the Rcpp package
library(Rcpp)

# Load the C++ implementation
sourceCpp("./code/cascade_estimation.cpp")


par_sweep_mf <- read_csv("./data/parsweep_noncc_4.csv")

no_sims <- 10000
net_2 <- newman_net(n = 5000, l = 0, m = 6)
adj_2 <- as_adj(net_2)

net_3 <- newman_net(n = 5000, l = 3, m = 0)
adj_3 <- as_adj(net_3)

net_4 <- four_cl_net_fn(n = 5000, m = 2)
adj_4 <- as_adj(net_4)


# sim_2 <- cascade_sim_par(1, net = net_2, adj = adj_2, p1 = 0.05, alpha = 0.1, total = no_sims)
# sim_2 <- sim_2 %>% mutate(net_size = 10000)
# 
# sim_3 <- cascade_sim_par(1, net = net_3, adj = adj_3, p1 = 0.05, alpha = 0.1, total = no_sims)
# sim_3 <- sim_3 %>% mutate(net_size = 10000)
# 
# sim_4 <- cascade_sim_par(1, net = net_4, adj = adj_4, p1 = 0.05, alpha = 0.1, total = no_sims)
# sim_4 <- sim_4 %>% mutate(net_size = 10000)

# compare -----------------------------------------------------------------

par_sweep_sim <- 
  bind_rows(
  tibble(p = seq(0.01, 0.19, length.out = 10), alpha = c(0)),
  tibble(p = seq(0.01, 0.24, length.out = 10), alpha = c(0.25)), 
  tibble(p = seq(0.01, 0.35, length.out = 10), alpha = c(0.5))
  ) %>% 
  mutate(res = list(NULL)) # %>% 
  # filter(
  #   (alpha == 0 & p <= 0.16) | 
  #     (alpha == 0.25 & p <= 0.22) | 
  #     (alpha == 0.5 & p <= 0.35)
  #   )


M <- 25000

for(i in 1:nrow(par_sweep_sim)){
  p <- par_sweep_sim$p[i]
  alpha <- par_sweep_sim$alpha[i]
  
  sim_2 <- run_simulations(adj_2 %>% as.matrix(), p, alpha, M) %>% 
    as_tibble() %>% mutate(network_type = "2")
  
  sim_3 <- run_simulations(adj_3 %>% as.matrix(), p, alpha, M) %>% 
    as_tibble() %>% mutate(network_type = "3")
  
  sim_4 <- run_simulations(adj_4 %>% as.matrix(), p, alpha, M) %>% 
    as_tibble() %>% mutate(network_type = "4")
  
  sum_df <- bind_rows(
    sim_2 %>% summarise(sim_es = mean(num_active), network_type = network_type[1]),
    sim_3 %>% summarise(sim_es = mean(num_active), network_type = network_type[1]),
    sim_4 %>% summarise(sim_es = mean(num_active), network_type = network_type[1])
  )
  
  par_sweep_sim$res[[i]] <- sum_df
  
  print(i/nrow(par_sweep_sim))
}



par_sweep_sim_sum <- par_sweep_sim %>% unnest(res) %>% # rename(network_type = networt_type) %>% 
  mutate(
    network_type = 
      case_when(network_type == "2" ~ "exp_size_2", 
                network_type == "3" ~ "exp_size_3", 
                network_type == "4" ~ "exp_size_4")
  )

test <- par_sweep_mf %>% filter(alpha %in% c(0,0.25, 0.5)) %>% 
  select(p, alpha, starts_with("exp")) %>% 
  pivot_longer(3:5, names_to = "network_type", values_to = "sim_es") %>% 
  filter(sim_es < 100)


test %>% filter(p <= 0.35) %>% 
  ggplot(aes(x = p, y = sim_es, color = factor(network_type)))  +
  geom_line() + 
  scale_y_log10() +
  geom_point(data = par_sweep_sim_sum) + 
  facet_wrap(~alpha, nrow = 3) 


test %>% filter(p <= 0.35) %>% 
  ggplot(aes(x = p, y = sim_es, color = factor(network_type)))  +
  geom_line() + 
  scale_y_log10() + 
  geom_point(data = par_sweep_sim_sum) + 
  facet_wrap(~alpha, nrow = 3) 


test %>% filter(p <= 0.35) %>% 
  ggplot(aes(x = p, y = sim_es, color = interaction(network_type,alpha)))  +
  geom_line() +
  scale_y_log10() +
  geom_point(data = par_sweep_sim_sum)




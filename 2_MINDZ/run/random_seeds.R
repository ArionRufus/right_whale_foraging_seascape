generate_seeds <- function() {
  # Generate random 31-bit parts
  seed1_part1 <- as.integer(sample(0:(2^31 - 1), 1))
  seed1_part2 <- as.integer(sample(0:(2^31 - 1), 1))
  seed2_part1 <- as.integer(sample(0:(2^31 - 1), 1))
  seed2_part2 <- as.integer(sample(0:(2^31 - 1), 1))
  
  # Convert the integers to hexadecimal format and combine them
  seed1_hex <- sprintf("Z'%08X%08X'", seed1_part1, seed1_part2)
  seed2_hex <- sprintf("Z'%08X%08X'", seed2_part1, seed2_part2)
  
  return(c(seed1_hex, seed2_hex))
}

# Generate and print seeds
seeds <- generate_seeds()
print(seeds)

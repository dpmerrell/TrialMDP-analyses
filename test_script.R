library("blockRARopt")

blockRARopt::block_rar_opt(64, 8, 0.5, 1.0, "results.sqlite", 0.5,0.5,0.5,0.5)# "beta_binom", "block_cost", "wald_failure")

print("About to connect")
conn = blockRARopt::connect_to_results("results.sqlite")

print("Successfully connected; about to fetch")
res = blockRARopt::fetch_result(conn, 0,0,0,0)

print("successfully fetched")
print(res)

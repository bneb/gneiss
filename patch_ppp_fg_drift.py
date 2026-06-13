with open("crates/gneiss-rtk/src/engine/ppp_fg.rs", "r") as f:
    content = f.read()

content = content.replace("let rcv_clk_drift = x_i[16];", "let rcv_clk_drift = if x_i.len() > 19 { x_i[19] } else { 0.0 };")
content = content.replace("x_i[0], x_i[1], x_i[2], x_i[3], x_i[4], x_i[5], x_i[16],", "x_i[0], x_i[1], x_i[2], x_i[3], x_i[4], x_i[5], if x_i.len() > 19 { x_i[19] } else { 0.0 },")

with open("crates/gneiss-rtk/src/engine/ppp_fg.rs", "w") as f:
    f.write(content)

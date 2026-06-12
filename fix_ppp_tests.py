import re

with open("crates/gneiss-rtk/src/engine/ppp.rs", "r") as f:
    content = f.read()

# remove anything outside mod tests if I messed up
# wait, earlier I ran: cat crates/gneiss-rtk/src/engine/ppp.rs | sed '/^#\[cfg(test)\]/,$d' > ppp.rs.tmp
# which REMOVED the tests. Then I appended scratch_test_ppp.rs, and then I appended MORE to the bottom.
# So I have two #[cfg(test)] blocks or one test block and some garbage at the end.

# Let's just find the first #[cfg(test)] and remove from there down.
idx = content.find("#[cfg(test)]")
if idx != -1:
    content = content[:idx]

with open("crates/gneiss-rtk/src/engine/ppp.rs", "w") as f:
    f.write(content)


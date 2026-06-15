with open("/Users/kevin/.gemini/antigravity/brain/5d333f93-5973-46d5-a672-726d6d8db288/task.md", "r") as f:
    text = f.read()

text = text.replace("- [ ] Write unit tests to cover new pure functions", "- [x] Write unit tests to cover new pure functions")
text = text.replace("- [ ] Eliminate mutants in updater_math.rs", "- [/] Eliminate mutants in updater_math.rs")

with open("/Users/kevin/.gemini/antigravity/brain/5d333f93-5973-46d5-a672-726d6d8db288/task.md", "w") as f:
    f.write(text)

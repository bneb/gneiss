with open("crates/gneiss-rtk/src/engine/measurement.rs", "r") as f:
    content = f.read()

replacement = """        state.add_ambiguity(ref_sat, 1, 5.0, 100.0);
        state.add_ambiguity(rov_sat1, 1, 10.0, 100.0);
        state.add_ambiguity(rov_sat2, 1, 15.0, 100.0);
        
        state.previous_windup.insert(ref_sat, 0.0);
        state.previous_windup.insert(rov_sat1, 0.0);
        state.previous_windup.insert(rov_sat2, 0.0);"""

content = content.replace("""        state.add_ambiguity(ref_sat, 1, 5.0, 100.0);
        state.add_ambiguity(rov_sat1, 1, 10.0, 100.0);
        state.add_ambiguity(rov_sat2, 1, 15.0, 100.0);""", replacement)

with open("crates/gneiss-rtk/src/engine/measurement.rs", "w") as f:
    f.write(content)

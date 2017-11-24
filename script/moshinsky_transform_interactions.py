
import mcscript

moshinsky_executable="/Users/amccoy/research/shell/programs/relutils/moshinsky"

def moshinsky_transform_operator():
	command_line=[
		moshinsky_executable,
		"< moshinsky.in"
	]

	mcscript.call(
        command_line,
        mode=mcscript.CallMode.kSerial
    )


def generate_moshinsky_tb_interactions(relative_interactions):
	for i in range(0,length(relative_interactions)):
		input_lines=[
			"tb 20 h2v15099",
			relative_interactions[i],
			"interaction{:04d}".format(i)
		]

		mcscript.utils.write_input("moshinsky.in",input_lines,verbose=True)
		moshinsky_transform_operator()


if (__name__ == "__MAIN__"):
	relative_interaction=[
		"jisp16_Nmax20_hw20.0_rel_Nmax10_u3st_1.0e-02.dat",
		"jisp16_Nmax20_hw20.0_rel_Nmax10_u3st_1.0e-03.dat",
		"jisp16_Nmax20_hw20.0_rel_Nmax10_u3st_1.0e-04.dat",
		"jisp16_Nmax20_hw20.0_rel_Nmax10_u3st_1.0e-05.dat",
		"jisp16_Nmax20_hw20.0_rel_Nmax20_N0_10.dat",
		"jisp16_Nmax20_hw20.0_rel_Nmax20_N0_12.dat",
		"jisp16_Nmax20_hw20.0_rel_Nmax20_N0_14.dat",
		"jisp16_Nmax20_hw20.0_rel_Nmax20_N0_16.dat",
		"jisp16_Nmax20_hw20.0_rel_Nmax20_N0_18.dat",
		"jisp16_Nmax20_hw20.0_rel_Nmax20_N0_20.dat",
		"jisp16_Nmax20_hw20.0_rel_Nmax20_Nrel_10.dat",
		"jisp16_Nmax20_hw20.0_rel_Nmax20_Nrel_12.dat",
		"jisp16_Nmax20_hw20.0_rel_Nmax20_Nrel_14.dat",
		"jisp16_Nmax20_hw20.0_rel_Nmax20_Nrel_16.dat",
		"jisp16_Nmax20_hw20.0_rel_Nmax20_Nrel_18.dat",
		"jisp16_Nmax20_hw20.0_rel_Nmax20_Nrel_20.dat",
		"jisp16_Nmax20_hw20.0_rel_Nmax20_u3st_1.0e-02.dat",
		"jisp16_Nmax20_hw20.0_rel_Nmax20_u3st_1.0e-03.dat",
		"jisp16_Nmax20_hw20.0_rel_Nmax20_u3st_1.0e-04.dat",
		"jisp16_Nmax20_hw20.0_rel_Nmax20_u3st_1.0e-05.dat"
	]

	generate_moshinsky_tb_interactions(relative_interactions)
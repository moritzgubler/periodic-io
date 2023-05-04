configure_file(test_all.sh test_all.sh)
add_subdirectory("test_dat")
add_test(NAME "Self_constistent_reading_writing" COMMAND "./test_all.sh"  "../src/transform/transform-periodic" "../src/fingerprint-programs/fingerprint-distance")
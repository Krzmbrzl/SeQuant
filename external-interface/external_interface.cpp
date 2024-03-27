#include "processing.hpp"

#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/parse_expr.hpp>
#include <SeQuant/core/utility/string.hpp>

#include <CLI/CLI.hpp>

#include <nlohmann/json.hpp>

#include <boost/algorithm/string.hpp>

#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <stdexcept>
#include <string>
#include <string_view>

using nlohmann::json;
using sequant::toUtf16;

void generate_itf(const json &blocks, std::string_view out_file) {
	for (const json &current_block : blocks) {
		const std::string block_name = current_block.at("name");

		for (const json &current_result : current_block.at("results")) {
			const std::string result_name = current_result.at("name");
			const std::string input_file  = current_result.at("equation_file");

			if (!std::filesystem::exists(input_file)) {
				throw std::runtime_error("Specified input file '" + input_file + "' does not exist");
			}

			// Read input file
			std::ifstream in(input_file);
			const std::string input(std::istreambuf_iterator< char >(in), {});

			const std::wstring transcoded_input = toUtf16(input);
			sequant::ExprPtr expression         = sequant::parse_expr(transcoded_input);

			expression = post_process(expression);

			// TODO: Add generated expression to ITF block to generate ITF code from
		}
	}
}

void generate_code(const json &details) {
	const std::string format   = details.at("output_format");
	const std::string out_path = details.at("output_path");

	if (boost::iequals(format, "itf")) {
		generate_itf(details.at("code_blocks"), out_path);
	} else {
		throw std::runtime_error("Unknown code generation target format '" + std::string(format) + "'");
	}
}

void process(const json &driver) {
	if (driver.contains("code_generation")) {
		const json &details = driver.at("code_generation");

		generate_code(details);
	}
}

int main(int argc, char **argv) {
	CLI::App app("Interface for reading in equations generated outside of SeQuant");
	argv = app.ensure_utf8(argv);

	std::string driver;
	app.add_option("--driver", driver, "Path to the JSON file used to drive the processing")->required();

	CLI11_PARSE(app, argc, argv);

	if (!std::filesystem::exists(driver)) {
		throw std::runtime_error("Specified driver file '" + driver + "' does not exist");
	}

	try {
		std::ifstream in(driver);
		json driver_info;
		in >> driver_info;

		process(driver_info);
	} catch (const std::exception &e) {
		std::wcout << "[ERROR]: " << e.what() << std::endl;
		return 1;
	}

	return 0;
}

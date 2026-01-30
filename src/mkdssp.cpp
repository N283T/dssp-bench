/*-
 * SPDX-License-Identifier: BSD-2-Clause
 *
 * Copyright (c) 2020 NKI/AVL, Netherlands Cancer Institute
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice, this
 *    list of conditions and the following disclaimer
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#if __has_include("config.hpp")
# include "config.hpp"
#endif

#include <chrono>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>

#include <cif++.hpp>
#include <mcfp/mcfp.hpp>
#include <nlohmann/json.hpp>

#include "dssp.hpp"
#include "revision.hpp"

namespace fs = std::filesystem;
using json = nlohmann::json;

// --------------------------------------------------------------------

// recursively print exception whats:
void print_what(const std::exception &e)
{
	std::cerr << e.what() << std::endl;
	try
	{
		std::rethrow_if_nested(e);
	}
	catch (const std::exception &nested)
	{
		std::cerr << " >> ";
		print_what(nested);
	}
}

// --------------------------------------------------------------------

// Check if the file is a JSON file based on extension
bool is_json_file(const fs::path &path)
{
	auto ext = path.extension();
	if (ext == ".gz" or ext == ".xz")
		ext = path.stem().extension();
	return ext == ".json";
}

// Load JSON and create a cif::datablock for DSSP processing
cif::file load_json_structure(const fs::path &path)
{
	// Read and decompress if needed
	std::string content;
	{
		cif::gzio::ifstream in(path);
		if (not in.is_open())
			throw std::runtime_error("Could not open JSON file");

		std::ostringstream ss;
		ss << in.rdbuf();
		content = ss.str();
	}

	auto j = json::parse(content);

	// Create a new CIF file and datablock
	cif::file f;
	auto &db = f.emplace_back("structure");

	// Create atom_site category
	auto &atom_site = db["atom_site"];

	// Create entity category (required for dssp)
	auto &entity = db["entity"];
	entity.emplace({
		{"id", "1"},
		{"type", "polymer"}
	});

	// Create struct_asym category (required for dssp)
	auto &struct_asym = db["struct_asym"];
	struct_asym.emplace({
		{"id", "A"},
		{"entity_id", "1"}
	});

	// Create entity_poly_seq for sequence information
	auto &entity_poly_seq = db["entity_poly_seq"];

	// Create pdbx_poly_seq_scheme for PDB numbering
	auto &pdbx_poly_seq_scheme = db["pdbx_poly_seq_scheme"];

	int atom_id = 1;
	const auto &residues = j["residues"];

	for (const auto &res : residues)
	{
		std::string chain_id = res["chain_id"].get<std::string>();
		int seq_id = res["seq_id"].get<int>();
		std::string compound_id = res["compound_id"].get<std::string>();

		// Add to entity_poly_seq
		entity_poly_seq.emplace({
			{"entity_id", "1"},
			{"num", std::to_string(seq_id)},
			{"mon_id", compound_id}
		});

		// Add to pdbx_poly_seq_scheme
		pdbx_poly_seq_scheme.emplace({
			{"asym_id", chain_id},
			{"entity_id", "1"},
			{"seq_id", std::to_string(seq_id)},
			{"mon_id", compound_id},
			{"pdb_strand_id", chain_id},
			{"pdb_seq_num", std::to_string(seq_id)},
			{"pdb_ins_code", "."}
		});

		// Add backbone atoms (N, CA, C, O) - required for DSSP
		const auto &atoms = res["atoms"];
		for (const auto &[atom_name, coords] : atoms.items())
		{
			double x = coords[0].get<double>();
			double y = coords[1].get<double>();
			double z = coords[2].get<double>();

			atom_site.emplace({
				{"id", std::to_string(atom_id++)},
				{"type_symbol", atom_name.substr(0, 1)},
				{"label_atom_id", atom_name},
				{"label_alt_id", "."},
				{"label_comp_id", compound_id},
				{"label_asym_id", chain_id},
				{"label_entity_id", "1"},
				{"label_seq_id", std::to_string(seq_id)},
				{"pdbx_PDB_ins_code", "."},
				{"Cartn_x", std::to_string(x)},
				{"Cartn_y", std::to_string(y)},
				{"Cartn_z", std::to_string(z)},
				{"occupancy", "1.00"},
				{"B_iso_or_equiv", "0.00"},
				{"pdbx_formal_charge", "."},
				{"auth_seq_id", std::to_string(seq_id)},
				{"auth_comp_id", compound_id},
				{"auth_asym_id", chain_id},
				{"auth_atom_id", atom_name},
				{"pdbx_PDB_model_num", "1"}
			});
		}

		// Add side chain atoms if present
		if (res.contains("side_chain"))
		{
			for (const auto &sc_atom : res["side_chain"])
			{
				std::string atom_name = sc_atom["name"].get<std::string>();
				const auto &pos = sc_atom["pos"];
				double x = pos[0].get<double>();
				double y = pos[1].get<double>();
				double z = pos[2].get<double>();

				// Determine element from atom name
				std::string element = atom_name.substr(0, 1);
				if (atom_name.size() > 1 && std::islower(atom_name[1]))
					element = atom_name.substr(0, 2);

				atom_site.emplace({
					{"id", std::to_string(atom_id++)},
					{"type_symbol", element},
					{"label_atom_id", atom_name},
					{"label_alt_id", "."},
					{"label_comp_id", compound_id},
					{"label_asym_id", chain_id},
					{"label_entity_id", "1"},
					{"label_seq_id", std::to_string(seq_id)},
					{"pdbx_PDB_ins_code", "."},
					{"Cartn_x", std::to_string(x)},
					{"Cartn_y", std::to_string(y)},
					{"Cartn_z", std::to_string(z)},
					{"occupancy", "1.00"},
					{"B_iso_or_equiv", "0.00"},
					{"pdbx_formal_charge", "."},
					{"auth_seq_id", std::to_string(seq_id)},
					{"auth_comp_id", compound_id},
					{"auth_asym_id", chain_id},
					{"auth_atom_id", atom_name},
					{"pdbx_PDB_model_num", "1"}
				});
			}
		}
	}

	// Add disulfide bonds if present
	if (j.contains("ss_bonds"))
	{
		auto &struct_conn = db["struct_conn"];
		int conn_id = 1;

		for (const auto &bond : j["ss_bonds"])
		{
			std::string chain1 = bond["chain1"].get<std::string>();
			int seq1 = bond["seq1"].get<int>();
			std::string chain2 = bond["chain2"].get<std::string>();
			int seq2 = bond["seq2"].get<int>();

			struct_conn.emplace({
				{"id", "disulf" + std::to_string(conn_id++)},
				{"conn_type_id", "disulf"},
				{"ptnr1_label_asym_id", chain1},
				{"ptnr1_label_comp_id", "CYS"},
				{"ptnr1_label_seq_id", std::to_string(seq1)},
				{"ptnr1_label_atom_id", "SG"},
				{"ptnr2_label_asym_id", chain2},
				{"ptnr2_label_comp_id", "CYS"},
				{"ptnr2_label_seq_id", std::to_string(seq2)},
				{"ptnr2_label_atom_id", "SG"},
				{"pdbx_ptnr1_PDB_ins_code", "."},
				{"pdbx_ptnr2_PDB_ins_code", "."}
			});
		}
	}

	return f;
}

// --------------------------------------------------------------------

int d_main(int argc, const char *argv[])
{
	using namespace std::literals;

	auto &config = mcfp::config::instance();

	config.init("Usage: mkdssp [options] input-file [output-file]",
		mcfp::make_option<std::string>("output-format", "Output format, can be either 'dssp' for classic DSSP or 'mmcif' for annotated mmCIF. The default is chosen based on the extension of the output file, if any."),
		mcfp::make_option<short>("min-pp-stretch", 3, "Minimal number of residues having PSI/PHI in range for a PP helix, default is 3"),
		mcfp::make_option("write-other", "If set, write the type OTHER for loops, default is to leave this out"),
		mcfp::make_option("no-dssp-categories", "If set, will suppress output of new DSSP output in mmCIF format"),

		mcfp::make_option("calculate-accessibility", "Default is to not calculate the surface accessibility when the output format is mmCIF"),

		mcfp::make_option<std::string>("mmcif-dictionary", "Path to the mmcif_pdbx.dic file to use instead of default"),

		mcfp::make_option("help,h", "Display help message"),
		mcfp::make_option("version", "Print version"),
		mcfp::make_option("verbose,v", "verbose output"),
		mcfp::make_option("quiet", "Reduce verbose output to a minimum"),
		mcfp::make_option("timing", "Show timing breakdown for benchmarking"),

		mcfp::make_hidden_option<int>("debug,d", "Debug level (for even more verbose output)"));

	config.parse(argc, argv);

	// --------------------------------------------------------------------

	if (config.has("version"))
	{
		write_version_string(std::cout, config.has("verbose"));
		exit(0);
	}

	if (config.has("help") or config.operands().empty())
	{
		std::cerr << config << std::endl;
		exit(config.has("help") ? 0 : 1);
	}

	if (config.has("output-format") and config.get<std::string>("output-format") != "dssp" and config.get<std::string>("output-format") != "mmcif")
	{
		std::cerr << "Output format should be one of 'dssp' or 'mmcif'" << std::endl;
		exit(1);
	}

	if (config.count("quiet"))
		cif::VERBOSE = -1;
	else
		cif::VERBOSE = config.count("verbose");

	// --------------------------------------------------------------------

	// private mmcif_pdbx dictionary?
	if (config.has("mmcif-dictionary"))
	{
		fs::path mmcif_dict = config.get<std::string>("mmcif-dictionary");

		cif::add_file_resource("mmcif_pdbx.dic", mmcif_dict);

		// Try to be smart, maybe dssp-extension.dic is at that location as well?
		auto dir = fs::canonical(mmcif_dict.parent_path());
		if (auto dssp_dict = cif::load_resource("dssp-extension.dic"); dssp_dict == nullptr and fs::exists(dir / "dssp-extension.dic"))
			cif::add_data_directory(dir);
	}

	cif::file f;
	fs::path input_path = config.operands().front();
	bool is_json = is_json_file(input_path);

	auto parse_start = std::chrono::high_resolution_clock::now();

	if (is_json)
	{
		// Load JSON format
		if (cif::VERBOSE > 0)
			std::cerr << "Loading JSON file...";

		f = load_json_structure(input_path);

		if (cif::VERBOSE > 0)
			std::cerr << " done\n";
	}
	else
	{
		// Load CIF/PDB format
		try
		{
			cif::gzio::ifstream in(input_path);
			if (not in.is_open())
			{
				std::cerr << "Could not open file" << std::endl;
				exit(1);
			}

			if (cif::VERBOSE > 0)
				std::cerr << "Loading file...";

			f.load(in);

			if (cif::VERBOSE > 0)
				std::cerr << " fixup file...";

			cif::pdb::fixup_pdbx(f);

			if (cif::VERBOSE > 0)
				std::cerr << " done\n";
		}
		catch (const std::exception &e)
		{
			std::cerr << e.what() << '\n';

			f = cif::pdb::read(input_path);
		}
	}

	auto parse_end = std::chrono::high_resolution_clock::now();

	// --------------------------------------------------------------------

	short pp_stretch = 3;
	if (config.has("min-pp-stretch"))
		pp_stretch = config.get<short>("min-pp-stretch");

	bool writeOther = config.has("write-other");

	std::string fmt;
	if (config.has("output-format"))
		fmt = config.get<std::string>("output-format");

	fs::path output;
	if (config.operands().size() > 1)
		output = config.operands()[1];

	if (fmt.empty() and not output.empty())
	{
		if (output.extension() == ".gz" or output.extension() == ".xz")
		{
			if (output.stem().extension() == ".dssp")
				fmt = "dssp";
			else
				fmt = "cif";
		}
		else if (output.extension() == ".dssp")
			fmt = "dssp";
		else
			fmt = "cif";
	}

	if (fmt == "dssp")
	{
		// See if the data will fit at all
		auto &db = f.front();
		for (const auto &[chain_id, seq_nr] : db["pdbx_poly_seq_scheme"].rows<std::string, int>("pdb_strand_id", "pdb_seq_num"))
		{
			if (chain_id.length() > 1 or seq_nr > 99999)
			{
				std::cerr << "The data in this file won't fit in the old DSSP format, please use the mmCIF format instead." << std::endl;
				exit(2);
			}
		}
	}

	auto dssp_start = std::chrono::high_resolution_clock::now();

	dssp dssp(f.front(), 1, pp_stretch, fmt == "dssp" or config.has("calculate-accessibility"));

	auto dssp_end = std::chrono::high_resolution_clock::now();

	if (config.has("timing"))
	{
		auto parse_ms = std::chrono::duration<double, std::milli>(parse_end - parse_start).count();
		auto dssp_ms = std::chrono::duration<double, std::milli>(dssp_end - dssp_start).count();
		auto stats = dssp.get_statistics();

		// Output in dssp-zig compatible format for benchmark runner
		std::cerr << "TIMING: residues=" << stats.count.residues
		          << " parse_ms=" << std::fixed << std::setprecision(3) << parse_ms
		          << " calc_total_ms=" << dssp_ms
		          << std::endl;
	}

	if (not output.empty())
	{
		cif::gzio::ofstream out(output);

		if (not out.is_open())
		{
			std::cerr << "Could not open output file" << std::endl;
			exit(1);
		}

		if (fmt == "dssp")
			dssp.write_legacy_output(out);
		else
		{
			dssp.annotate(f.front(), writeOther, not config.has("no-dssp-categories"));
			out << f.front();
		}
	}
	else
	{
		if (fmt == "dssp")
			dssp.write_legacy_output(std::cout);
		else
		{
			dssp.annotate(f.front(), writeOther, not config.has("no-dssp-categories"));
			std::cout << f.front();
		}
	}

	return 0;
}

// --------------------------------------------------------------------

int main(int argc, const char *argv[])
{
	int result = 0;

	try
	{
#if defined(DATA_DIR)
		cif::add_data_directory(DATA_DIR);
#endif
		result = d_main(argc, argv);
	}
	catch (const std::exception &ex)
	{
		print_what(ex);
		exit(1);
	}

	return result;
}

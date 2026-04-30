use clap::Parser;
use cli::Subcommands;

mod bam_concordance;
mod cli;
mod common;
fn main() {
    let args = cli::Cli::parse();

    match args.command {
        Subcommands::BamConcordance(bam_concordance_args) => {
            bam_concordance::bam_concordance(&bam_concordance_args).unwrap();
        }
        _ => panic!("Not Implemented yet"),
    }
}

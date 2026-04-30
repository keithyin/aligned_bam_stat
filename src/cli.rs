
use clap::{self, Parser, Subcommand, Args};

#[derive(Parser)]
pub struct Cli {

    #[command(subcommand)]
    pub command: Subcommands

}


#[derive(Subcommand)]
pub enum Subcommands {
    BamConcordance(BamConcordanceArgs)
}


#[derive(Args)]
pub struct BamConcordanceArgs {
    pub reffasta: String,
    pub aligned_bam: String,

    #[arg(long="hcregions")]
    pub hcregions: Option<String>,
    #[arg(long="hcvariants")]
    pub hcvariants: Option<String>,
    
    pub chrom: Option<String>,

}
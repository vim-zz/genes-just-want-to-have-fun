// START
// Generate the initial population
// Compute fitness
// REPEAT
//     Selection
//     Crossover
//     Mutation
//     Compute fitness
// UNTIL population has converged
// STOP

use rand::{
    distributions::{Distribution, Uniform, Bernoulli},
    seq::IteratorRandom,
    Rng,
};

use std::fmt;

fn main() {
    let target = "hello yehonatan my son!";

    let chromosome_size = target.len();
    let general_population_size = 1000;
    let breeders_pairs = 10;
    let mutation_factor = 0.1;
    let fitness_factor = 1.0;

    let fitness = (chromosome_size as f64 * -fitness_factor) as isize;
    println!("fitness: {}", fitness);
    let mut rng = rand::thread_rng();
    let mut population = vec![];
    for _ in 0..general_population_size {
        population.push(Chromosome::from("_______________________"));
    }
    let mut alltimes_best_score = 0;

    for generation in 0.. {
        // Selection
        population.sort_by_key(|x| x.score(target));
        let breeders: Vec<_> = population.iter().take(breeders_pairs * 2).collect();

        let scores: Vec<_>= population.iter().map(|x| x.score(target)).collect();
        let avg_score = scores.iter().sum::<isize>() / population.len() as isize;
        let best_score = scores.iter().min().unwrap().to_owned();
        if best_score < alltimes_best_score || generation % 20 == 0 {
            let best = breeders.iter().find(|x| x.score(target) == best_score).unwrap();
            println!("=============================================================================");
            println!("best: {}, generation: {} ....... avg_score: {} ....... best_score: {}", best, generation, avg_score, best_score);
            println!("=============================================================================");
            for i in &breeders {
                println!("{}", i);
            }

            alltimes_best_score = best_score;
        }

        let crossover_from = rng.gen_range(0, chromosome_size);
        let crossover_to = rng.gen_range(crossover_from, chromosome_size);
        let crossover = (crossover_from, crossover_to);
        let (boys, girls): (Vec<_>, Vec<_>) = breeders
            .chunks_exact(2)
            .map(|x| breed(&x[0], &x[1], crossover))
            .unzip();

        // Breeding
        let mut children: Vec<_> = boys.into_iter()
            .chain(girls.into_iter())
            .collect();

        children.iter_mut().for_each(|x| x.mutate(mutation_factor));

        // Survival
        let weakest = general_population_size - breeders_pairs * 2;
        let _: Vec<_>= population.drain(weakest..).collect();
        population.append(&mut children);

        let qualified: Vec<_> = population.iter().filter(|x| x.score(target) < fitness).collect();
        if qualified.len() > 0 {
            println!("DONE after {} genrations", generation);
            println!("Qualified population: {:?}", qualified);
            break;
        }
    }

}

fn breed(papa: &Chromosome, mama: &Chromosome, crossover: (usize, usize)) -> (Chromosome, Chromosome) {
    // crossover
    let (from, to) = crossover;
    let mixed_genes: (Vec<Gene>, Vec<Gene>) = papa.genes.iter()
        .zip(mama.genes.iter())
        .enumerate()
        .map(|(i, (a, b))| {
            if i >= from && i <= to {
                (b, a)
            } else {
                (a, b)
            }
        })
        .unzip();

    let boy = Chromosome{genes: mixed_genes.0};
    let girl = Chromosome{genes: mixed_genes.1};
    println!("mama+papa: {} + {} => {} + {}", mama, papa, boy, girl);
    (boy, girl)
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
struct Gene(char);

impl Gene {
    fn new(ch: char) -> Self {
        Gene(ch)
    }

    fn random<R: Rng + ?Sized>(rng: &mut R) -> Gene {
        let ch = rng.sample_iter(SourceCode).next().unwrap();
        Gene(ch)
    }

    fn mutate<R: Rng + ?Sized>(&mut self, rng: &mut R, p: f64) {
        let bernoulli = Bernoulli::new(p).unwrap();
        if bernoulli.sample(rng) {
            *self = Gene::random(rng);
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Hash, PartialOrd, Ord)]
struct Chromosome {
    pub genes: Vec<Gene>,
}

impl Chromosome {
    fn new(genes: Vec<Gene>) -> Self {
        Chromosome {
            genes,
        }
    }

    fn random(len: usize) -> Self {
        let mut rng = rand::thread_rng();
        let mut genes = vec![];
        for _ in 0..len {
            genes.push(Gene::random(&mut rng))
        }
        Chromosome { genes }
    }

    fn score<T: Into<Chromosome>>(&self, target: T) -> isize {
        let sum: isize = self.genes 
            .iter()
            .zip(target.into().genes.iter())
            .map(|(a, b)| if a == b { 1 } else { 0 })
            .sum();

        sum * -1
    }

    fn mutate(&mut self, p: f64) {
        let mut rng = rand::thread_rng();
        self.genes.iter_mut().choose(&mut rng).unwrap().mutate(&mut rng, p);
    }
}

impl From<&str> for Chromosome {
    fn from(txt: &str) -> Chromosome {
        Chromosome {
            genes: txt.chars().map(|x| Gene(x)).collect()
        }
    }
}

impl fmt::Display for Chromosome {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.genes.iter().map(|x| x.0).collect::<String>())
    }
}

struct SourceCode;

impl Distribution<char> for SourceCode {
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> char {
        const GEN_SRCCODE_CHARSET: &[u8] = b"ABCDEFGHIJKLMNOPQRSTUVWXYZ\
                abcdefghijklmnopqrstuvwxyz\
                0123456789\
                +-*/%^!&|<>@_.,;:#$?(){}[] ";
        let range = GEN_SRCCODE_CHARSET.len();
        let range = Uniform::new(0, range);
        let i = rng.sample_iter(range).next().unwrap();
        return GEN_SRCCODE_CHARSET[i as usize] as char;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test] 
    fn mutate_chromosome() {
        let mut c9s1 = Chromosome::from("____________________");
        let c9s2 = c9s1.clone();
        c9s1.mutate(1.0);
        println!("original: {}, mutation: {}", c9s1, c9s2);
        assert_ne!(c9s1, c9s2);
    }

    #[test]
    fn chromosomes_are_eq_by_value() {
        assert_eq!(Chromosome::from("x"), Chromosome::from("x"));
    }

    #[test]
    fn random_chromosomes() {
        let chromosome_size = 5;
        let mut chromosomes = vec![];
        for _ in 0..3 {
            chromosomes.push(Chromosome::random(chromosome_size));
        }
        println!("random chromosomes: {:#?}", chromosomes);
        assert_ne!(chromosomes[0], chromosomes[1]);
        assert_ne!(chromosomes[1], chromosomes[2]);
    }

    #[test]
    fn sorted_chromosomes() {
        let mut chromosomes: Vec<Chromosome> = vec![];

        chromosomes.push("hxxxx".into());
        chromosomes.push("hellx".into());
        chromosomes.push("hexxx".into());

        // Selection
        println!("before: {:#?}", chromosomes);
        chromosomes.sort_by_key(|x| x.score("hello"));
        println!("after: {:#?}", chromosomes);

        assert_eq!(chromosomes[0].score("hello"), -5);
        assert_eq!(chromosomes[2].score("hello"), 0);
    }
}

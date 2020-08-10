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

fn main() {
    let mut rng = rand::thread_rng();

    let chromosome_size = 500;
    let fitness = (chromosome_size as f64 * -0.9) as isize;
    println!("fitness: {}", fitness);

    let general_population_size = 1000;
    let mut population = vec![Chromosome::new(chromosome_size); general_population_size];
    let mut alltimes_best_score = 0;

    for genration in 0..10000 {
        // Selection
        population.sort_by_key(|x| x.score());
        let breeders_pairs = 2;
        let breeders: Vec<_> = population.iter().take(breeders_pairs * 2).collect();

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

        children.iter_mut().for_each(|x| x.mutate());

        // Survival
        let weakest = general_population_size - breeders_pairs * 2;
        let _: Vec<_>= population.drain(weakest..).collect();
        population.append(&mut children);

        let scores: Vec<_>= population.iter().map(|x| x.score()).collect();
        let avg_score = scores.iter().sum::<isize>() / population.len() as isize;
        let best_score = scores.iter().min().unwrap().to_owned();
        if best_score < alltimes_best_score {
            println!("generation: {} ....... avg_score: {} ....... best_score: {}", genration, avg_score, best_score);
            alltimes_best_score = best_score;
        }

        let qualified: Vec<_> = population.iter().filter(|x| x.score() < fitness).collect();
        if qualified.len() > 0 {
            println!("DONE after {} genrations", genration);
            println!("Qualified population: {:?}", qualified);
            break;
        }
    }

}

fn breed(papa: &Chromosome, mama: &Chromosome, crossover: (usize, usize)) -> (Chromosome, Chromosome) {
    // crossover
    // println!("mama {:?}", mama);
    // println!("papa {:?}", papa);
    let (from, to) = crossover;
    let (boy, girl): (Vec<Gene>, Vec<Gene>) = papa.genes.iter()
        .zip(mama.genes.iter())
        // .inspect(|x| println!("about to breed: {:?}", x))
        .enumerate()
        .map(|(i, (a, b))| {
            if i >= from && i <= to {
                (b, a)
            } else {
                (a, b)
            }
        })
        // .inspect(|x| println!("mazal tov! {:?}", x))
        .unzip();

    (Chromosome{genes: boy}, Chromosome{genes: girl})
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
enum Gene {
    Zero,
    One,
}

impl Gene {
    fn random<R: Rng + ?Sized>(rng: &mut R) -> Gene {
        let uniform = Uniform::new_inclusive(0, 1);
        match uniform.sample(rng) {
            0 => Gene::Zero,
            _ => Gene::One,
        }
    }

    fn mutate<R: Rng + ?Sized>(&mut self, rng: &mut R) {
        let bernoulli = Bernoulli::new(0.1).unwrap();
        if bernoulli.sample(rng) {
            *self = match self {
                Gene::Zero => Gene::One,
                Gene::One => Gene::Zero,
            };
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Hash, PartialOrd, Ord)]
struct Chromosome {
    pub genes: Vec<Gene>,
}

impl Chromosome {
    fn new(len: usize) -> Self {
        Chromosome {
            genes: vec![Gene::Zero; len],
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

    fn score(&self) -> isize {
        let sum: isize = self
            .genes
            .iter()
            .map(|x| match x {
                Gene::Zero => 0,
                Gene::One => 1,
            })
            .sum();

        sum * -1
    }

    fn mutate(&mut self) {
        let mut rng = rand::thread_rng();
        self.genes.iter_mut().choose(&mut rng).unwrap().mutate(&mut rng);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test] 
    fn mutate_chromosome() {
        let mut c9s = Chromosome::new(100);
        assert_eq!(c9s.score(), 0);
        c9s.mutate();
        println!("mutation: {:?}", c9s);
        assert_ne!(c9s.score(), 0);
    }

    #[test]
    fn chromosomes_are_eq_by_value() {
        let chromosome_size = 5;
        let mut c9s = Chromosome::new(chromosome_size);
        c9s.genes = vec![Gene::Zero; chromosome_size];
        assert_eq!(Chromosome::new(chromosome_size), c9s);
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
        let chromosome_size = 5;
        let mut chromosomes = vec![];

        let mut c9s = Chromosome::new(chromosome_size);
        c9s.genes = vec![Gene::Zero, Gene::Zero, Gene::Zero, Gene::Zero, Gene::Zero];
        chromosomes.push(c9s);

        let mut c9s = Chromosome::new(chromosome_size);
        c9s.genes = vec![Gene::One, Gene::One, Gene::One, Gene::One, Gene::One];
        chromosomes.push(c9s);

        let mut c9s = Chromosome::new(chromosome_size);
        c9s.genes = vec![Gene::Zero, Gene::Zero, Gene::Zero, Gene::Zero, Gene::Zero];
        chromosomes.push(c9s);

        // Selection
        println!("before: {:#?}", chromosomes);
        chromosomes.sort_by_key(|x| x.score());
        println!("after: {:#?}", chromosomes);

        assert_eq!(chromosomes[0].score(), -5);
        assert_eq!(chromosomes[2].score(), 0);
    }
}

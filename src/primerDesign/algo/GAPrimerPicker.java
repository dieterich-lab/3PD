/**
 * 
 */
package primerDesign.algo;

import java.util.Random;

import primerDesign.dsc.PrimerPair;
import primerDesign.dsc.RestrictionSite;
import primerDesign.util.Constants;
import primerDesign.util.SimpleTimer;
import cern.colt.list.ObjectArrayList;

/**
 * This class implements a genetic algorithm to pick the best primer set.
 * THIS IS A DEADEND CLASS - NOT TO BE FURTHER IMPLEMENTED!!!
 * 
 * This methodology is based on: Wu et.al.: "Primer design for multiplex PCR using a genetic algorithm. Soft Comp (2007) 11:855-863
 * 
 * @author Sebastian Fršhler
 *
 */
public class GAPrimerPicker{ //implements PrimerPickingAlgorithm{
	private ObjectArrayList chromosomes;
	private GAChromosome currentChromosome;
	private PrimerPair currentPrimerPair;
	private Random random;
	private double bestWeakness;
	private double weaknessSum;
	private boolean printDebugOutput = false;
	private SimpleTimer timer = new SimpleTimer();
	
	/**
	 * This method implements a primer picking strategy using a genetic algorithm.
	 * 
	 * @param optimalSites an array of optimal restriction sites, each associated with valid candidate primers to build a best set from
	 */
	public PrimerPair[] pickBestPrimerSet(RestrictionSite[] optimalSites) {		
		random = new Random();
		int randomUpstreamPrimerIndex;
		int randomDownstreamPrimerIndex;
		int randomTaqManProbeIndex;
		weaknessSum = 0;
		this.bestWeakness = Double.MAX_VALUE;
		this.chromosomes = new ObjectArrayList();
	
		if(printDebugOutput) System.out.print("Initializing population");
		// generate initial population 'P' by random sampling
		for(int i=0; i<Constants.GA_NUM_INDIVIDUALS; i++){
			currentChromosome = new GAChromosome(optimalSites);
			for(int j=0; j<optimalSites.length; j++){
				currentPrimerPair = new PrimerPair();
				
				randomUpstreamPrimerIndex = random.nextInt(optimalSites[j].getNumberOfValidUpstreamPrimers());
				randomDownstreamPrimerIndex = random.nextInt(optimalSites[j].getNumberOfValidDownstreamPrimers());
				randomTaqManProbeIndex = random.nextInt(optimalSites[j].getNumberOfValidTaqManProbes());
				
				currentPrimerPair.setForwardPrimer(optimalSites[j].getUpstreamPrimer(randomUpstreamPrimerIndex));
				currentPrimerPair.setReversePrimer(optimalSites[j].getDownstreamPrimer(randomDownstreamPrimerIndex));
				currentPrimerPair.setHybridizationProbe(optimalSites[j].getTaqManProbe(randomTaqManProbeIndex));
				
				currentChromosome.addTargetRegion(currentPrimerPair);
			}
			this.chromosomes.add(currentChromosome);
		}
		this.computeChromosomeWeaknessValues();
		
		if(printDebugOutput) System.out.println(" - done in " + timer.getTimeString());
		
		// set mutation and crossover probability -> defined in 'Constants'
		
		if(printDebugOutput) System.out.print("Picking best primer set");
		// while termination condition is NOT satified
		int iterations = Constants.GA_NUM_GENERATIONS;
		double currentLeastWeakness = 0;
		int chromosomes = this.chromosomes.size();
		
		while(iterations > 0 && this.bestWeakness - currentLeastWeakness >= Constants.GA_EPSILON){
			currentLeastWeakness = Double.MAX_VALUE;
		
			// for each element in the current population
			for(int i=0; i< chromosomes; i++){
			
				// randomly select two chromosomes (=candidate primer sets) from population
				GAChromosome first = this.rouletteWheel(random);
				GAChromosome second = this.rouletteWheel(random);
				
				// clone those two chromosomes - primers and sites are NOT cloned!!! since only references are switched during mutation and crossover
				// only PrimerSets need to be cloned
				GAChromosome new1 = first.shallowClone();
				GAChromosome new2 = second.shallowClone();

				// mutate each of the two cloned chromosomes
				new1.mutate();
				new2.mutate();
		
				// perform crossover of the two chromosomes
				new2 = new1.crossover(new2);
		
				// add the cloned cromosomes into mating pool
				this.chromosomes.add(new1);
				this.chromosomes.add(new2);
			}
			
			// select |P| top-scoring chromosomes to replace the original population (= to enter next round XOR be returned)
			// compute weakness of each chromosome
			this.computeChromosomeWeaknessValues();
			
			// order chromosomes by ascending weakness
			this.chromosomes.quickSort();
			
			// remember best weakness value
			currentLeastWeakness = ((GAChromosome)this.chromosomes.getQuick(0)).weakness;
			if(currentLeastWeakness < this.bestWeakness) this.bestWeakness = currentLeastWeakness;
			
			// retain |P| best chromosomes
			this.chromosomes = this.chromosomes.partFromTo(0, Constants.GA_NUM_INDIVIDUALS -1);
			
			iterations--;
		}
		
		if(printDebugOutput) System.out.println(" - done in " + timer.getTimeString());
		
		GAChromosome bestChromosome = (GAChromosome)this.chromosomes.getQuick(0);
		PrimerPair[] bestPairs = new PrimerPair[bestChromosome.targetRegions.size()];
		for(int i=0; i<bestChromosome.targetRegions.size(); i++){
			bestPairs[i] = (PrimerPair)bestChromosome.targetRegions.getQuick(i);
		}
			
		return bestPairs;
	}
	
	private void computeChromosomeWeaknessValues(){
		this.weaknessSum = 0;
		for(int i=0; i<this.chromosomes.size(); i++){
			this.weaknessSum += ((GAChromosome)this.chromosomes.getQuick(i)).calculateWeakness();
		}
	}
	
	private GAChromosome rouletteWheel(Random random){
		double sum = 0;
		double value = random.nextDouble();
		GAChromosome currentChromosome = null;
		
		for (int i=0; i<this.chromosomes.size(); i++) {
			currentChromosome = (GAChromosome)this.chromosomes.getQuick(i);
			if(sum/weaknessSum < value && (sum += currentChromosome.weakness)/weaknessSum >= value) break;
		}
		return currentChromosome;
	}
	
	private class GAChromosome implements Comparable{
		// vector containing primer pairs for each targte region
		private ObjectArrayList targetRegions;
		private int changePrimerType;
		private int changePosition;
		private boolean changeLongArm;
		private RestrictionSite[] sites;
		private double weakness;
		private double epsilonWeakness = 1E-03;
		private boolean hasBeenModified;

		GAChromosome(RestrictionSite[] sites) {
			this.sites = sites;
			this.targetRegions = new ObjectArrayList();
			this.hasBeenModified = true;
		}
		
		GAChromosome(RestrictionSite[] sites, ObjectArrayList targetRegions) {
			this.sites = sites;
			this.targetRegions = targetRegions;
			this.hasBeenModified = true;
		}
		
		void addTargetRegion(PrimerPair pair){
			this.targetRegions.add(pair);
			this.hasBeenModified = true;
		}
		
		void mutate(){
			// random element from same position and same primer type
			changePrimerType = random.nextInt(3);
			changePosition = random.nextInt(targetRegions.size());
			
			if(printDebugOutput) System.out.println("Mutating chromosome:\nBefore:\n " + this.toString());
			
			assert(changePosition < this.targetRegions.size());
			
			switch(changePrimerType){ 
				case 0 : ((PrimerPair)this.targetRegions.getQuick(changePosition)).setForwardPrimer(this.sites[changePosition].getUpstreamPrimer(random.nextInt(this.sites[changePosition].getNumberOfValidUpstreamPrimers()))) ; break;
				case 1 : ((PrimerPair)this.targetRegions.getQuick(changePosition)).setReversePrimer(this.sites[changePosition].getDownstreamPrimer(random.nextInt(this.sites[changePosition].getNumberOfValidDownstreamPrimers()))) ; break;
				case 2 : ((PrimerPair)this.targetRegions.getQuick(changePosition)).setHybridizationProbe(this.sites[changePosition].getTaqManProbe(random.nextInt(this.sites[changePosition].getNumberOfValidTaqManProbes()))) ; break;
				default : throw new IllegalStateException("Unhandled case!");
			}
			this.hasBeenModified = true;
			
			if(printDebugOutput) System.out.println("After:\n" + this.toString());
		}
		
		GAChromosome crossover(GAChromosome other){
			changePosition = random.nextInt(targetRegions.size());
			changeLongArm = random.nextBoolean();
			Object temp;
			
			assert(changePosition < this.targetRegions.size());
			
			if(printDebugOutput) System.out.println("Crossing over chromosome at position " + changePosition + ":\nBefore:\n " + this.toString() + "\n" + other.toString());
			
			if(changeLongArm){
				for(int i=changePosition; i>=0; i--){
					temp = this.targetRegions.getQuick(i);
					this.targetRegions.setQuick(i, other.targetRegions.getQuick(i));
					other.targetRegions.setQuick(i, temp);
				}
			}else{
				for(int i=changePosition; i< targetRegions.size(); i++){
					temp = this.targetRegions.getQuick(i);
					this.targetRegions.setQuick(i, other.targetRegions.getQuick(i));
					other.targetRegions.setQuick(i, temp);
				}
			}
			this.hasBeenModified = true;
			other.hasBeenModified = true;
			
			if(printDebugOutput) System.out.println("After:\n" + this.toString() + "\n" + other.toString());
			
			return other;
		}
		
		double calculateWeakness(){
			
			double weakness = 0;
			
			for(int i=0; i<this.targetRegions.size()-1; i++){
				for(int j=i+1; j<this.targetRegions.size(); j++){
					weakness += ((PrimerPair)this.targetRegions.getQuick(i)).getAlignmentValues(((PrimerPair)this.targetRegions.getQuick(j)));
				}
			}
			this.weakness = weakness;
			return this.weakness;
		}
		
		double getWeakness(){
			if(!hasBeenModified){
				return weakness;
			}
			else{
				return this.calculateWeakness();
			}
		}
		
		ObjectArrayList copyTargetRegions(){
			ObjectArrayList result = new ObjectArrayList();
			
			for(int i=0; i<this.targetRegions.size();i++){
				result.add(new PrimerPair((PrimerPair)this.targetRegions.getQuick(i)));
			}
			
			return result;
		}
		
		GAChromosome shallowClone(){
			return new GAChromosome(this.sites, this.copyTargetRegions());
		}

		public int compareTo(Object o) {
			GAChromosome other = (GAChromosome) o;
			if(this.weakness < other.weakness - this.epsilonWeakness) return -1;
			else if(this.weakness > other.weakness + this.epsilonWeakness) return 1;
			else return 0;
		}
		
		public String toString(){
			StringBuffer buffer = new StringBuffer();
			for(int i=0; i<this.targetRegions.size(); i++){
				buffer.append(((PrimerPair)this.targetRegions.getQuick(i)).toString() + "\n");
			}
			return buffer.toString();
		}
	}
}

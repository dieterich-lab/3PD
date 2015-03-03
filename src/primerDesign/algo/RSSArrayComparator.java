package primerDesign.algo;

import java.util.Comparator;

import primerDesign.dsc.RestrictionSite;

public class RSSArrayComparator implements Comparator {

	public int compare(Object o1, Object o2) {
		RestrictionSite first = (RestrictionSite) o1;
		RestrictionSite second = (RestrictionSite) o2;
		if(first.getNumberOfValidUpstreamPrimers() < second.getNumberOfValidUpstreamPrimers()) return -1;
		else if(first.getNumberOfValidUpstreamPrimers() > second.getNumberOfValidUpstreamPrimers()) return 1;
		else if(first.getNumberOfValidUpstreamPrimers() == second.getNumberOfValidUpstreamPrimers()) return 0;
		else throw new IllegalStateException("Unhandled case");
	}

}

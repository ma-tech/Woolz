package uk.ac.mrc.hgu.TVImage;

public class NodeInfo
{
	String nodeID;
	String emageID;
	String correlation;

	public NodeInfo(String nodeID)
	{
		this.nodeID = nodeID;
	}

	public static boolean isLeaf(String id)
	{
		return id.startsWith("GENE");
	}

	public static String getEmageId(String id)
	{
		return null;
	}

	public static String getURLForGeneSearch(String id)
	{
		String url = "http:// genex.hgu.mrc.ac.uk/das/jsp/geneByComponent.jsp?tissue=EMAP:" + id;
		return url;
	}

	public static String getURLForGeneEmageId(String id)
	{
		String url = "http:// genex.hgu.mrc.ac.uk/das/jsp/geneByComponent.jsp?tissue=EMAP:" + id;
		return url;
	}
}

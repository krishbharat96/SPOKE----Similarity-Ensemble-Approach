from neo4j.v1 import GraphDatabase, basic_auth
import parse_sea as ps
import json

CYPHER_Add_Assay_Edge = """
    MATCH (p:Protein), (c:Compound) WHERE p.chembl_id = '{p_id}'
    AND c.chembl_id = '{c_id}'
    CREATE (c)-[r:INTERACTS_CiP]->(p) SET r.assays = {assay_arr},
    r.source = 'ChEMBL', r.action_type = 'UNKNOWN', r.license = 'CC BY-SA 3.0',
    r.drug_interaction = False
"""

CYPHER_Add_Assay_No_Edge = """
    MATCH (p:Protein)<-[r:INTERACTS_CiP]-(c:Compound) WHERE
    p.chembl_id = '{p_id}' AND c.chembl_id = '{c_id}'
    SET r.assays = {assay_arr}
"""

CYPHER_Update_Assay_Edge = """
    MATCH (p:Protein)<-[r:INTERACTS_CiP]-(c:Compound)
    WHERE p.chembl_id = '{p_id}' AND c.chembl_id = '{c_id}'
    SET r.assays = {assay_arr}
"""

CYPHER_Vestige_Assay_Edge = """
    MATCH (c:Compound)-[r:INTERACTS_CiP]->(p:Protein)
    WHERE c.chembl_id = '{c_id}' AND p.chembl_id = '{p_id}'
    SET r.assays = null
"""

CYPHER_Remove_Excess_Assays = """
    MATCH (c:Compound)-[r:INTERACTS_CiP]->(p:Protein)
    WHERE not exists (r.assays) AND r.drug_interaction = False
    AND exists (c.chembl_id) AND exists (p.chembl_id)
    AND r.action_type = 'UNKNOWN' DELETE r
"""

CYPHER_Obtain_Assays = """
    MATCH (p:Protein)<-[r:INTERACTS_CiP]-(c:Compound)
    WHERE c.identifier contains '' AND p.identifier contains ''
    AND exists (r.assays)
    RETURN p.chembl_id as p_id, c.chembl_id as c_id,
    r.assays as r_assays
"""

CYPHER_Obtain_All_Edges = """
    MATCH (p:Protein)<-[r:INTERACTS_CiP]-(c:Compound)
    WHERE c.identifier contains '' AND p.identifier contains ''
    AND exists(c.chembl_id) AND exists(p.chembl_id)
    RETURN p.chembl_id as p_id, c.chembl_id as c_id
"""

def main():
    print "Connecting to GraphDB..."
    driver = GraphDatabase.driver("bolt://127.0.0.1/:7687", auth=basic_auth("neo4j", "neo4j2"))
    session = driver.session()
    print "Connected!"
    print "Obtaining old assay values..."
    old_recs = session.run(CYPHER_Obtain_Assays)
    old_dict = dict()
    for o_r in old_recs:
        concat_id = str(o_r["c_id"]) + "|" + str(o_r["p_id"])
        if not concat_id in old_dict.keys():
            old_dict.update({concat_id:o_r["r_assays"]})
    print "Obtained!"
    print "Parsing new assay file"
    new_dict = ps.get_assay_dict_chembl('full_assays_act.tsv')
    print "Parsed!"
    
    print "Updating Assays.."
    update_assays(new_dict, old_dict, session)
    session.close()
           
def update_assays(new_dict, old_dict, session):
    count_vestige = 0
    count_add = 0
    count_update = 0

    all_edge_arr = []
    all_edges = session.run(CYPHER_Obtain_All_Edges)

    for ae in all_edges:
        concat_id = str(ae["c_id"]) + "|" + str(ae["p_id"])
        if not concat_id in all_edge_arr:
            all_edge_arr.append(concat_id)
    
    vestige_arr = list(set(old_dict.keys()) - set(new_dict.keys()))
    add_arr = list(set(new_dict.keys()) - set(old_dict.keys()))

    for v in vestige_arr:
        cmpd_id = v.split("|")[0]
        prot_id = v.split("|")[1]
        session.run(CYPHER_Vestige_Assay_Edge.format(c_id=cmpd_id,
                                                     p_id=prot_id))
        count_vestige = count_vestige + 1
        old_dict.pop(v)

    for a in add_arr:
        cmpd_id = a.split("|")[0]
        prot_id = a.split("|")[1]
        if not a in all_edge_arr:
            session.run(CYPHER_Add_Assay_Edge.format(c_id=cmpd_id,
                                                     p_id=prot_id,
                                                     assay_arr=str(new_dict[a])))
            all_edge_arr.append(a)
        else:
            session.run(CYPHER_Add_Assay_No_Edge.format(c_id=cmpd_id,
                                                        p_id=prot_id,
                                                        assay_arr=str(new_dict[a])))
            
        count_add = count_add + 1
        new_dict.pop(a)

    for u in new_dict.keys():
        cmpd_id = u.split("|")[0]
        prot_id = u.split("|")[1]
        if not new_dict[u] == old_dict[u]:
            session.run(CYPHER_Update_Assay_Edge.format(c_id=cmpd_id,
                                                        p_id=prot_id,
                                                        assay_arr=str(new_dict[u])))
            count_update = count_update + 1

    # session.run(CYPHER_Remove_Excess_Assays)

    print "COUNT UPDATED : " + str(count_update)
    print "COUNT ADDED : " + str(count_add)
    print "COUNT VESTIGES : " + str(count_vestige)

main()






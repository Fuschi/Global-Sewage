# Tables Directory Overview

This directory contains summary tables that provide insights and statistical data from the network analysis of ARG and mOTU.

### Table: info_nodes
**Filename**: `info_nodes.tsv`  
**Description**: This table contains detailed information about the nodes in the network, including their community memberships and specific metadata related to either mOTU or AMR genes. The table is organized to distinguish between microbial taxa and antimicrobial resistance genes effectively.

**Columns**:
- **mgnet**: Indicates the specific network analysis the node is part of (e.g., mOTU+resfinder or mOTU+functional).
- **taxa_id**: General identifier for a node, which can be an mOTU or an AMR gene.
- **comm_id**: The community to which the node belongs.
- **source**: Identifies whether the node is a bacterium or an AMR gene.
- **type**: Specifies the node type, such as bacterium, resfinder, or functional gene.

**mOTU Metadata** (applicable only to bacterial nodes, `NA` for AMR nodes):
- **path_taxs**
- **kingdom**
- **phylum**
- **class**
- **order**
- **family**
- **genus**
- **species_gtdb**
- **human_gut**
- **is_from_gut**

**AMR Metadata** (applicable only to AMR nodes, `NA` for mOTU nodes):
- **group_amr**
- **biocide_class_amr**
- **class_amr**
- **cluster_representative**
- **cluster_version**
- **metal_class_amr**
- **resistance_type**

### Table 2: Edges Data Frame
**Filename**: `info_edges_all.tsv`  
**Description**: This table encapsulates the relationships between nodes in the network, highlighting the interactions between microbial taxa and antimicrobial resistance genes. It provides a comprehensive view of the connectivity and the community crossing events within the network.

**Columns**:
- **mgnet**: Identifies the network analysis (e.g., mOTU+resfinder or mOTU+functional) from which the edge data is derived.
- **weight**: Quantitative measure of the connection correlation of the CLR abundances between two nodes.
- **node_i**, **node_j**: Identifiers for the nodes at either end of an edge.
- **link_type**: Describes the nature of the link (e.g. mOTU-ResFinder, mOTU-mOTU, mOTU-Functional...).
- **crossing_communities**: Indicates whether the edge crosses between different microbial communities.

**Metadata for Node i** (suffix `_i` indicates data for `node_i`):
- **source_i**: Whether the node is a bacterium or an AMR gene.
- **type_i**: Type of node (bacterium, resfinder, functional gene).
- **path_taxs_i**, **kingdom_i**, **phylum_i**, **class_i**, **order_i**, **family_i**, **genus_i**, **species_gtdb_i**: Taxonomic metadata for mOTU nodes.
- **human_gut_i**, **is_from_gut_i**: Indicators of human gut origin.
- **group_amr_i**, **biocide_class_amr_i**, **class_amr_i**, **cluster_representative_i**, **cluster_version_i**, **metal_class_amr_i**, **resistance_type_i**: Metadata related to AMR characteristics.

**Metadata for Node j** (suffix `_j` indicates data for `node_j`).

### Table 3: info_edges_only_mOTU_AMR
**Filename**: `info_edges_only_mOTU_AMR.tsv`  
**Description**: This table is a specialized subset of the `edges` data frame. It focuses specifically on the interactions between microbial taxa (mOTU) and antimicrobial resistance genes (ARG), capturing only the relationships where the correlations are positive.

\name{gsea.rename.genesets.edb}
\alias{gsea.rename.genesets.edb}
\title{Parse an edb XML object, and rename the genesets with optional prefix and...}
\usage{gsea.rename.genesets.edb(edb, old.names, new.names)
}
\description{Parse an edb XML object, and rename the genesets with optional prefix and
suffix.}
\details{This is intended to be used by gsea.rename.genesets, which should calculate
old.names and new.names. Entries in edb, given that they are xml are in no
strict order (although they seem to be ordered from most -ve NES to most
+ve NES). Thus you need to supply the old.names is so that we can ensure
the entries in edb get renamed correctly.}
\value{An mutable XMLInternalDOM (ie an XML tree much like an XMLDocument
but which can be modified) which is a copy of the original edb
XMLDocument.}
\author{Mark Cowley, 2009-10-06}
\arguments{\item{edb}{an XMLDocument. See import.gsea.edb(..., edb=TRUE)}
\item{old.names}{a vector of original geneset names. these should be found
within the edb names}
\item{new.names}{the new names, in the same order as the old names.}
}


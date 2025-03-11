"""
Useful bits of code that do not belong in anywhere else and should be imported
"""

def add_annotation(soup, tag, key, value):
    """
    Add as an annotation the key value pair to the given tag.
    """
    # Find list of annotations or add one.
    annotation = tag.find("annotation")
    if annotation is None:
        annotation = soup.new_tag('annotation')
        tag.append(annotation)
    rdf = annotation.find("rdf:RDF")
    if rdf is None:
        rdf = soup.new_tag('rdf:RDF', attrs={
            "xmlns:bqbiol":"http://biomodels.net/biology-qualifiers/",
            "xmlns:bqmodel":"http://biomodels.net/model-qualifiers/",
            "xmlns:dcterms":"http://purl.org/dc/terms/",
            "xmlns:rdf":"http://www.w3.org/1999/02/22-rdf-syntax-ns#",
            "xmlns:vCard":"http://www.w3.org/2001/vcard-rdf/3.0#",
            "xmlns:vCard4":"http://www.w3.org/2006/vcard/ns#"})
        annotation.append(rdf)
    description = rdf.find('rdf:Description')
    if description is None:
        description = soup.new_tag('rdf:Description',
            attrs={'rdf:about':f'#{tag.get("metaid","NONE")}'})
        rdf.append(description)
    bqbiol = description.find('bqbiol:is')
    if bqbiol is None:
        bqbiol = soup.new_tag('bqbiol:is')
        description.append(bqbiol)
    holder = bqbiol.find('rdf:Bag')
    if holder is None:
        holder = soup.new_tag('rdf:Bag')
        bqbiol.append(holder)

    # Add the annotation in an rdf part which is 'holder'
    holder.append(soup.new_tag('rdf:li',attrs={
        'rdf:resource': f'https://identifiers.org/{key}/{value}'
    }))

def main():
    """
    A busy way to do nothing.
    """

if __name__ == '__main__':
    main()

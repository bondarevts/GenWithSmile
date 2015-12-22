from __future__ import absolute_import
import pytest

from gws.io.star_smiles import parse_star_smiles
from gws.io.star_smiles.star_smiles_parser import ParserResults


@pytest.mark.parametrize('star_smiles,expected_result', [
    ('NCC(=O)O', ParserResults(
        smiles='NCC(=O)O',
        attach_bonds=[],
        attach_positions=[],
        insert_positions=[],
        fragment_pos=[]
    )),
    ('NCl*C', ParserResults(
        smiles='NClC',
        attach_bonds=[],
        attach_positions=[2],
        insert_positions=[],
        fragment_pos=[]
    )),
    ('NCl**C', ParserResults(
        smiles='NClC',
        attach_bonds=[],
        attach_positions=[],
        insert_positions=[2],
        fragment_pos=[]
    )),
    ('NCl***C', ParserResults(
        smiles='NClC',
        attach_bonds=[],
        attach_positions=[2],
        insert_positions=[2],
        fragment_pos=[]
    )),
    ('NCl<-=#>***C', ParserResults(
        smiles='NClC',
        attach_bonds=[(2, '-'), (2, '='), (2, '#')],
        attach_positions=[2],
        insert_positions=[2],
        fragment_pos=[]
    )),
    ('NCl^C', ParserResults(
        smiles='NClC',
        attach_bonds=[],
        attach_positions=[],
        insert_positions=[],
        fragment_pos=[2]
    )),
    ('{NCl}^C', ParserResults(
        smiles='NClC',
        attach_bonds=[],
        attach_positions=[],
        insert_positions=[],
        fragment_pos=[0, 1, 2]
    )),
    ('{NCl}C', ParserResults(
        smiles='NClC',
        attach_bonds=[],
        attach_positions=[0, 1, 2],
        insert_positions=[0, 1, 2],
        fragment_pos=[]
    )),
    ('C{N}{Cl}C', ParserResults(
        smiles='CNClC',
        attach_bonds=[],
        attach_positions=[1, 2, 3],
        insert_positions=[1, 2, 3],
        fragment_pos=[]
    )),
    ('C{N}*{Cl}**C', ParserResults(
        smiles='CNClC',
        attach_bonds=[],
        attach_positions=[1],
        insert_positions=[2, 3],
        fragment_pos=[]
    )),
    ('{CC}^CC^CC**<#>C<=>***C{CC<->C}{C}*C<->{CC}**C*', ParserResults(
        smiles='CCCCCCCCCCCCCCCC',
        attach_bonds=[(5, '#'), (6, '='), (9, '-'), (12, '-')],
        attach_positions=[6, 8, 9, 10, 11, 15],
        insert_positions=[5, 6, 8, 9, 10, 13, 14],
        fragment_pos=[0, 1, 3]
    ))
])
def test_without_marks(star_smiles, expected_result):
    result = parse_star_smiles(star_smiles)
    assert result == expected_result

echo $#
if (($# == 0))
then
   echo "1 : maillage"
   echo "2 : image fonctionnelle"
   echo "3 : proj sphere 5"
   echo "4 : proj sphere 3"
   echo "5 : proj ligne ext"
   echo "6 : fichier final"
else
   echo "# tree 1.0

*BEGIN TREE EXECUTE

*BEGIN TREE CreateControlWindow

*END

*END" > aux1
(more aux1; echo "# tree 1.0

*BEGIN TREE EXECUTE

*BEGIN TREE LoadObject") > aux2
(more aux2; echo "filename    $1") > aux1
(more aux1; echo "res_pointer 1

*END

*END") > aux2
(more aux2; echo "# tree 1.0

*BEGIN TREE EXECUTE

*BEGIN TREE LoadObject") > aux1
(more aux1; echo "filename    $2") > aux2
(more aux2; echo "res_pointer 2

*END

*END")>aux1
(more aux1; echo "# tree 1.0

*BEGIN TREE EXECUTE

*BEGIN TREE FusionObjects
objects     2 1
res_pointer 3
method      Fusion3DMethod

*END

*END")>aux2
(more aux2; echo "# tree 1.0

*BEGIN TREE EXECUTE

*BEGIN TREE FusionObjects
objects     2 1
res_pointer 4
method      Fusion3DMethod

*END

*END")>aux1
(more aux1; echo "# tree 1.0

*BEGIN TREE EXECUTE

*BEGIN TREE FusionObjects
objects     2 1
res_pointer 5
method      Fusion3DMethod

*END

*END") > aux2
(more aux2; echo "# tree 1.0

*BEGIN TREE EXECUTE

*BEGIN TREE Fusion3DParams
object      3
method      sphere
sub_method  mean
depth       5
step        2.5

*END

*END")>aux1
(more aux1; echo "# tree 1.0

*BEGIN TREE EXECUTE

*BEGIN TREE Fusion3DParams
object      4
method      sphere
sub_method  mean
depth       3
step        1.5

*END

*END")> aux2
(more aux2 ; echo "# tree 1.0

*BEGIN TREE EXECUTE

*BEGIN TREE Fusion3DParams
object      5
method      line_external
sub_method  mean
depth       3
step        1.5

*END

*END") > aux1
(more aux1; echo "# tree 1.0

*BEGIN TREE EXECUTE

*BEGIN TREE ExportTexture")>aux2
(more aux2;echo "filename $3") > aux1
(more aux1;echo "object   3

*END

*END") >aux2
(more aux2; echo "# tree 1.0

*BEGIN TREE EXECUTE

*BEGIN TREE ExportTexture")>aux1
(more aux1;echo "filename $4") > aux2
(more aux2;echo "object   4

*END

*END")>aux1
(more aux1; echo "# tree 1.0

*BEGIN TREE EXECUTE

*BEGIN TREE ExportTexture")>aux2
(more aux2;echo "filename $5")>aux1
(more aux1;echo "object   5

*END

*END") > $6
rm -f aux1
rm -f aux2

fi
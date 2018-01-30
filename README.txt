build.sh permet de compiler les differentes versions sans fftw.

build_with_fftw.sh compile les differentes versions avec fftw, avec la verification des valeurs activee.
une variable fftw_path permet de preciser ou trouver fftw si elle est compilee localement.

Chaque binaire peut etre appele avec en parametre le fichier wav a lire et le nombre de samples a lire.
measure.py peut etre appele avec le binaire a mesurer. Il faut qu'un fichier long_audio.wav d'au moins
2^27 samples soit present dans le repertoire. Un tel fichier est telechargeable ici :
https://bigfiles.univ-paris8.fr/2bzvxgd3

Seuls les wav mono 16bit unsigned little endian sont supportes (pas tous).

Le fichier callgrind.out est le dump de profilage mentionne dans le rapport.
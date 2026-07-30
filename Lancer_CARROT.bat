@echo off
REM ═══════════════════════════════════════════════════════════════════════════
REM Lancer CARROT.bat — Lanceur automatique pour Windows
REM ═══════════════════════════════════════════════════════════════════════════
REM Double-cliquez sur ce fichier pour lancer CARROT avec un acces complet a
REM votre dossier utilisateur (%USERPROFILE%), sans taper la moindre commande.
REM
REM Ce que fait ce script :
REM   - Monte votre dossier utilisateur (%USERPROFILE% : Bureau, Documents,
REM     Telechargements...) dans le conteneur, sous /data : shinyFiles y verra
REM     tous vos vrais dossiers, et AutoSpectral pourra y lire vos fichiers FCS
REM     et y ecrire ses resultats directement, sans copie intermediaire.
REM   - Supprime automatiquement tout ancien conteneur "carrot" avant d'en
REM     relancer un nouveau, pour eviter une erreur si vous relancez ce
REM     script plusieurs fois.
REM   - Ouvre automatiquement votre navigateur une fois l'application prete.
REM
REM Prerequis : Docker Desktop installe et lance (icone de la baleine active
REM dans la barre des taches). Le partage de fichiers de Docker Desktop gere
REM automatiquement les permissions sur Windows : pas besoin de --user ici.
REM ═══════════════════════════════════════════════════════════════════════════

set IMAGE=camellialambert/carrot
set NOM_CONTENEUR=carrot
set PORT_HOTE=80

echo ================================================
echo   Lancement de CARROT
echo ================================================

REM Verifie que Docker repond
docker version >nul 2>&1
if errorlevel 1 (
    echo Docker ne repond pas. Verifiez que Docker Desktop est bien lance
    echo ^(icone de la baleine dans la barre des taches^), puis reessayez.
    echo Telechargement : https://docs.docker.com/get-docker/
    pause
    exit /b 1
)

REM Supprime un eventuel ancien conteneur CARROT
docker rm -f %NOM_CONTENEUR% >nul 2>&1

echo Verification de la derniere version de CARROT...
docker pull %IMAGE%

echo Demarrage du conteneur...
docker run -dp %PORT_HOTE%:3838 --rm --name %NOM_CONTENEUR% -v "%USERPROFILE%":/data %IMAGE%

echo Attente du demarrage de l'application...
timeout /t 8 /nobreak >nul

set URL=http://localhost:%PORT_HOTE%/
echo CARROT est pret : %URL%

start "" "%URL%"

echo.
echo Pour arreter CARROT : docker stop %NOM_CONTENEUR%
pause

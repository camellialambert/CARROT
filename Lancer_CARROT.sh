#!/bin/bash
# ═══════════════════════════════════════════════════════════════════════════
# Lancer CARROT.sh — Lanceur automatique pour Linux et macOS
# ═══════════════════════════════════════════════════════════════════════════
# Double-cliquez sur ce fichier (ou exécutez-le dans un terminal) pour lancer
# CARROT avec un accès complet à votre dossier personnel ($HOME), sans avoir
# à taper la moindre commande Docker vous-même.
#
# Ce que fait ce script :
#   - Monte votre dossier personnel ($HOME, ex: Bureau, Documents...) dans le
#     conteneur, sous /data : shinyFiles y verra tous vos vrais dossiers, et
#     AutoSpectral pourra y lire vos fichiers FCS et y écrire ses résultats
#     directement, sans copie intermédiaire.
#   - Supprime automatiquement tout ancien conteneur "carrot" avant d'en
#     relancer un nouveau, pour éviter une erreur "name already in use" si
#     vous relancez ce script plusieurs fois.
#   - Ouvre automatiquement votre navigateur une fois l'application prête.
#
# NOTE — --platform linux/amd64 : NÉCESSAIRE sur Mac Apple Silicon (M1/M2/M3),
# car l'image a été construite pour amd64 (x86_64). Docker Desktop l'exécute
# alors via une émulation Rosetta interne. Sans impact sur Linux/Mac Intel,
# où le mode natif est déjà amd64.
#
# NOTE — PAS de --user ici (contrairement à une version précédente de ce
# script) : l'image rocker/shiny a besoin de démarrer en tant que root en
# interne (via s6-overlay) pour ensuite lancer Shiny Server correctement.
# Imposer un utilisateur arbitraire dès le départ empêche ce démarrage et
# fait planter le conteneur en boucle (erreur vue précédemment : "The user
# 'undefined' does not have permissions..."). Docker Desktop (Mac/Windows)
# gère l'écriture dans le dossier monté normalement malgré l'absence de
# --user, sans configuration supplémentaire.
# ═══════════════════════════════════════════════════════════════════════════

set -e

IMAGE="camellialambert/carrot"
NOM_CONTENEUR="carrot"
PORT_HOTE=80

echo "════════════════════════════════════════════"
echo "  Lancement de CARROT"
echo "════════════════════════════════════════════"

# Vérifie que Docker est installé et accessible
if ! command -v docker &> /dev/null; then
    echo "❌ Docker n'est pas installé ou n'est pas dans le PATH."
    echo "   Installez Docker Desktop : https://docs.docker.com/get-docker/"
    read -p "Appuyez sur Entrée pour fermer..."
    exit 1
fi

# Supprime un éventuel ancien conteneur CARROT (pour éviter un conflit de nom)
if docker ps -a --format '{{.Names}}' | grep -q "^${NOM_CONTENEUR}$"; then
    echo "→ Arrêt du conteneur CARROT précédent…"
    docker rm -f "$NOM_CONTENEUR" > /dev/null 2>&1 || true
fi

# Récupère la dernière version de l'image (ne fait rien si déjà à jour)
echo "→ Vérification de la dernière version de CARROT…"
docker pull --platform linux/amd64 "$IMAGE"

# Lance le conteneur :
#   --platform linux/amd64 -> force l'architecture de l'image (voir note ci-dessus)
#   -v "$HOME":/data        -> expose tout votre dossier personnel dans le conteneur
#   --rm                    -> nettoie le conteneur automatiquement à l'arrêt
echo "→ Démarrage du conteneur…"
docker run --platform linux/amd64 -dp ${PORT_HOTE}:3838 \
  --rm \
  --name "$NOM_CONTENEUR" \
  -v "$HOME":/data \
  "$IMAGE"

# Attend que Shiny Server réponde avant d'ouvrir le navigateur
echo "→ Attente du démarrage de l'application…"
for i in $(seq 1 30); do
    if curl -s -o /dev/null "http://localhost:${PORT_HOTE}/"; then
        break
    fi
    sleep 1
done

URL="http://localhost:${PORT_HOTE}/"
echo "✔ CARROT est prêt : $URL"

# Ouvre le navigateur automatiquement (macOS: open, Linux: xdg-open)
if command -v open &> /dev/null; then
    open "$URL"          # macOS
elif command -v xdg-open &> /dev/null; then
    xdg-open "$URL"      # Linux
else
    echo "Ouvrez manuellement votre navigateur à l'adresse : $URL"
fi

echo ""
echo "Pour arrêter CARROT : docker stop $NOM_CONTENEUR"
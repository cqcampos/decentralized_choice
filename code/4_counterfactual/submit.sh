#!/bin/bash
# Description: Wrapper to submit Slurm jobs with account and user-specific email
set -euo pipefail

# --- Configuration ---
SEND_EMAIL=false   # override at runtime: SEND_EMAIL=false bash submit.sh

# --- Context ---
ACCOUNT=${SLURM_JOB_ACCOUNT:-pi-cqcampos}
USER_NAME=${USER:-$(whoami)}

# --- Map usernames/accounts to email addresses ---
case "$ACCOUNT:$USER_NAME" in
  faculty:*)      EMAIL_USER="Christopher.Campos@chicagobooth.edu" ;;
  *:atumturk)     EMAIL_USER="Ayse.Tumturk@chicagobooth.edu" ;;
  *:cfogal)       EMAIL_USER="Connor.Fogal@chicagobooth.edu" ;;
  *:ryanlee22)    EMAIL_USER="Ryan.Lee2@chicagobooth.edu" ;;
  *)              EMAIL_USER="${USER_NAME}@chicagobooth.edu" ;;
esac

# Make path for logs in the current working directory if it doesn't already exist
mkdir -p logs

# --- Run the main scripts ---
echo "Submitting job:"
echo "  Account: $ACCOUNT"
echo "  User: $USER_NAME"

if [ "$SEND_EMAIL" = true ]; then
  echo "  Email notifications: enabled ($EMAIL_USER)"
    sbatch --account="$ACCOUNT" \
        --mail-user="$EMAIL_USER" \
        --mail-type=END,FAIL \
        <(sed '1s/^\xEF\xBB\xBF//; s/\r$//' counterfactual.sh)

else
  echo "  Email notifications: disabled"
      sbatch --account="$ACCOUNT" \
        <(sed '1s/^\xEF\xBB\xBF//; s/\r$//' counterfactual.sh)
fi
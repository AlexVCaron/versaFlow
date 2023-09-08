#!/usr/bin/env bash


# If Tower token (TOWER_ACCESS_TOKEN) and workspace ID (TOWER_WORKSPACE_ID) are set, fill the config file
if [ -n "$TOWER_ACCESS_TOKEN" ] && [ -n "$TOWER_WORKSPACE_ID" ]; then

echo "Tower token and workspace ID are set, setting config file"
cat <<EOF > .base/tower.config
tower {
    enabled = true
    accessToken = '$TOWER_ACCESS_TOKEN'
    workspaceId = '$TOWER_WORKSPACE_ID'
}
EOF

fi

package com.averycorp.prismtask.data.remote.sync

import kotlinx.coroutines.flow.Flow
import kotlinx.coroutines.flow.MutableStateFlow
import kotlinx.coroutines.flow.StateFlow
import kotlinx.coroutines.flow.asStateFlow
import kotlinx.coroutines.flow.map
import kotlinx.coroutines.flow.update

enum class SyncState {
    IDLE,
    PENDING,
    IN_PROGRESS,
    SUCCESS,
    ERROR,
}

data class SyncStatus(
    val state: SyncState = SyncState.IDLE,
    val lastSyncTime: Long = 0L,
    val pendingChanges: Int = 0,
    val errorMessage: String? = null,
)

class SyncStateRepository {

    private val _syncStatus = MutableStateFlow(SyncStatus())

    // StateFlow already guarantees distinct emissions — no distinctUntilChanged needed.
    val syncStatus: StateFlow<SyncStatus> = _syncStatus.asStateFlow()

    val isSyncing: Flow<Boolean> = syncStatus.map { it.state == SyncState.IN_PROGRESS }

    val hasPendingChanges: Flow<Boolean> = syncStatus.map { it.pendingChanges > 0 }

    fun markSyncStarted() {
        _syncStatus.update { it.copy(state = SyncState.IN_PROGRESS, errorMessage = null) }
    }

    fun markSyncSuccess(syncTime: Long = System.currentTimeMillis()) {
        _syncStatus.update {
            it.copy(
                state = SyncState.SUCCESS,
                lastSyncTime = syncTime,
                pendingChanges = 0,
                errorMessage = null,
            )
        }
    }

    fun markSyncError(message: String) {
        _syncStatus.update { it.copy(state = SyncState.ERROR, errorMessage = message) }
    }

    fun markSyncPending(pendingCount: Int = 1) {
        _syncStatus.update {
            it.copy(state = SyncState.PENDING, pendingChanges = pendingCount)
        }
    }

    fun incrementPendingChanges() {
        _syncStatus.update { it.copy(pendingChanges = it.pendingChanges + 1) }
    }

    fun decrementPendingChanges() {
        _syncStatus.update {
            val remaining = (it.pendingChanges - 1).coerceAtLeast(0)
            it.copy(
                pendingChanges = remaining,
                state = if (remaining == 0 && it.state == SyncState.PENDING) SyncState.IDLE
                        else it.state,
            )
        }
    }

    fun reset() {
        _syncStatus.value = SyncStatus()
    }
}

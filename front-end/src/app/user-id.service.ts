import { BehaviorSubject, Observable } from 'rxjs';
import { Injectable } from '@angular/core';

import { BackendHttpClient } from './backend-http.client';
import { UserIdServiceInterface } from './user-id.service.interface';

@Injectable({
  providedIn: 'root'
})
export class UserIdService implements UserIdServiceInterface {

  private readonly userId$: BehaviorSubject<string | null> = new BehaviorSubject<string | null> (null);
  readonly userId: Observable<string | null> = this.userId$.asObservable();

  constructor(private backendHttpClient: BackendHttpClient) { 
    this.setUserId();
  }

  /**
   * Asks the backendHttpClient to get the user id, and updates the user id when it receives a response.
   */
  async setUserId(): Promise<void> {
    const userId = await this.backendHttpClient.getUserId(localStorage.getItem('userId') || '');
    localStorage.setItem('userId', userId);
    this.userId$.next(userId);
  }
}

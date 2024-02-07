import { BehaviorSubject, Observable } from 'rxjs';
import { Injectable } from '@angular/core';

import { UserIdServiceInterface } from './user-id.service.interface';

@Injectable({
  providedIn: 'root'
})
export class MockUserIdService implements UserIdServiceInterface {

  private readonly userId$: BehaviorSubject<string | null> = new BehaviorSubject<string | null> (null);
  readonly userId: Observable<string | null> = this.userId$.asObservable();

  async setUserId(): Promise<void> {
    this.userId$.next('mock_user_id');
  }
}
